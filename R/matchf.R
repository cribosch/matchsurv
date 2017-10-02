### {{{ compdata
##' Data structured for matchpropexc
##' @param formula formula with 'Surv' outcome (see \code{coxph}); avoid using strata() or factor()
##' @param data data frame
##' @param idControl vector control indicator (idControl==1 indicates exposed individual in cluster i)
##' @param cluster vector cluster indicator (one cluster for each exposed individual)
##' @examples 
##' dd<-data.sim(nca=5000, ncont=5)
##' setdd<-compdata(Surv(time, status)~x+z+cc, cluster=id, idControl=j, data=dd)
##' summary(setdd)
##' @author Cristina Boschini
##' @return A setup dataset, ready for \code{matchpropexc}
##' @export
compdata<-function(formula, data, cluster, idControl,...){
  currentOPTs <- options("na.action")
  options(na.action = "na.pass")
  cl <- match.call()
  m <- match.call(expand.dots=TRUE)[1:5]
  Terms <- terms(formula,data=data, idControl=idControl, cluster=cluster)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m,"response")
  if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
  if (ncol(Y)==2) {
    exit <- Y[,1]
    entry <- rep(0,nrow(Y))
    status <- Y[,2]
    Truncation <- FALSE
  } else {
    entry <- Y[,1]
    exit <- Y[,2]
    status <- Y[,3]
    Truncation <- TRUE
  }
  # cluster <- NULL
  # if (!is.null(attributes(Terms)$specials$cluster)){
  #   ts <- survival::untangle.specials(Terms, "cluster")
  #   Terms <- Terms[-ts$terms]
  #   cluster <- m[[ts$vars]]
  # }
  idControl<-model.extract(m,"idControl")
  names(idControl)<- NULL
  cluster<-model.extract(m,"cluster")
  names(cluster)<- NULL
  
  X <- model.matrix(Terms, m)
  options(na.action = currentOPTs$na.action)
  
  if (!is.null(intpos <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  p<-ncol(X)
  if(!is.null(colnames(X))) namesX<-colnames(X)
  else if(p>0) namesX <- paste("var",seq(1,p),sep="")
  
  if (Truncation){
    
    wp <- mets::familyclusterWithProbands.index(cluster, idControl, Rindex==1)
    
    start<- apply(cbind(entry[wp$pairs[,1]],entry[wp$pairs[,2]]),1,max)
    end<- apply(cbind(exit[wp$pairs[,1]],exit[wp$pairs[,2]]),1,min)
    indicator <- ifelse(end==exit[wp$pairs[,2]],
                        status[wp$pairs[,2]],ifelse(status[wp$pairs[,1]]>0,-1,0))
    group <- cluster[wp$pairs[,1]]
    subj <- idControl[wp$pairs[,1]]-1
    size.group <- mets::cluster.index(group)$cluster.size
    
    funw <- function(x) sum(x==1)
    
    weight <- rep(tapply(indicator,group,funw), size.group)
    weight <- c(ifelse(indicator==-1,-1,ifelse(indicator==0,1,weight)))
    
    fund <- function(x) (!(duplicated(x>1) & x>1))*x
    
    weight <- unlist(tapply(weight,group,fund))
    stat <- ifelse(indicator==-1,1,
                   ifelse(indicator==1,
                          ifelse(weight>=1,1,0),0))
    
    d3 <- data.frame(start,end,stat,group,subj,weight)
    colnames(d3)<- c("entry","exit", "status","cluster","unexp.subj","weight")
    
    if (ncol(X)>0) {
      Xcases <- as.matrix(X[wp$pairs[,2],], ncol=ncol(X))
      colnames(Xcases)<-colnames(X)
      d3 <- data.frame(d3,Xcases, check.names=FALSE)
    }
    d3<-d3[d3$start<=d3$end,]
  } else {
    
    wp <- mets::familyclusterWithProbands.index(cluster,idControl,Rindex=1)
    end<- apply(cbind(exit[wp$pairs[,1]],exit[wp$pairs[,2]]),1,min)
    indicator <- ifelse(end==exit[wp$pairs[,2]],
                        status[wp$pairs[,2]],ifelse(status[wp$pairs[,1]]>0,-1,0))
    group <- cluster[wp$pairs[,1]]
    subj <- idControl[wp$pairs[,1]]-1
    size.group <- mets::cluster.index(group)$cluster.size
    
    funw <- function(x) sum(x==1)
    
    weight <- rep(tapply(indicator,group,funw), size.group)
    weight <- c(ifelse(indicator==-1,-1,ifelse(indicator==0,1,weight)))
    
    fund <- function(x) (!(duplicated(x>1) & x>1))*x
    
    weight <- unlist(tapply(weight,group,fund))
    stat <- ifelse(indicator==-1,1,
                   ifelse(indicator==1,
                          ifelse(weight>=1,1,0),0))
    
    d3 <- data.frame(end,stat,group,subj,weight)
    colnames(d3)<- c("exit","status","cluster","unexp.subj","weight")
    
    if (ncol(X)>0) {
      Xcases <- as.matrix(X[wp$pairs[,2],], ncol=ncol(X))
      colnames(Xcases)<-colnames(X)
      d3 <- data.frame(d3,Xcases, check.names=FALSE)
    }
  }
  return(d3)
}


###}}} compdata

###{{{ matchpropexc0
matchpropexc0 <- function(X,entry, exit, status, weight,
                          strata=NULL, beta,stderr=TRUE,
                          strata.name=NULL,...){
  if(is.vector(X)) X <- matrix(X, ncol=1)
  p <-ncol(X)
  if (missing(beta)) beta <-rep(0,p)
  if (p==0) X<-cbind(rep(0,length(exit)))
  if (!is.null(strata)){
    nstrata<-nlevels(strata)
    stratalev <- levels(strata)
    strataidx <- lapply(stratalev, function(x) which(strata==x))
    if (!all(unlist(lapply(strataidx, function(x) length(x)>0))))
      stop("Strata without any observation")
    dd <- lapply(strataidx, function(ii)
      .Call("prep",
            entry[ii], exit[ii], status[ii], weight[ii],
            as.matrix(X)[ii,,drop=FALSE],sum(entry[ii])!=0,
            package="matchsurv"))
    
    obj <- function(pp, U=FALSE, all=FALSE) {
      val <- lapply(dd, function (d)
        with(d,
             .Call("PL",pp,X,XX,Sign,jumps,weight, package="matchsurv")))
      gradient <- Reduce("+", lapply(val, function(x) x$gradient))
      hessian <- Reduce("+", lapply(val, function(x) x$hessian))
      S0 <- lapply(val, function(x) x$S0)
      nevent<-unlist(lapply(S0,length))
      if (all){
        U <- do.call("rbind", lapply(val, function(x) x$U))
        time<-lapply(dd, function(x) x$time[x$ord+1])
        ord<-lapply(dd, function(x) x$ord+1)
        jumps<-lapply(dd, function(x)  x$jumps+1)
        jumpstime<-lapply(dd, function(x) x$time[x$ord+1][x$jumps+1])
        weight<-lapply(val, function(x)  x$weight)
        E<-lapply(val, function(x) x$E)
        xjumps<-lapply(val, function(x) x$xjumps)
        #S0 <- lapply(val, function(x) x$S0)
        #nevent<-unlist(lapply(S0,length))
        return(list(gradient=gradient, hessian=hessian,
                    U=U, S0=S0, nevent=nevent,
                    ord=ord, time=time, jumps=jumps, 
                    jumpstime=jumpstime, weight=weight,
                    E=E, xjumps=xjumps))
      }
      structure(nevent,gradient=-gradient, hessian=-hessian)
    }
  } else {
    nstrata<-1
    dd <- .Call("prep",
                entry,exit, status,weight,X,
                sum(entry)!=0,
                package="matchsurv")
    
    obj <- function(pp, U=FALSE, all=FALSE) {
      val <- with(dd,
                  .Call("PL",pp,X,XX,Sign,jumps,weight, package="matchsurv"))
      val$nevent<-length(val$S0)
      if (all){
        val$time<-dd$time[dd$ord+1]
        val$ord<-dd$ord+1
        val$jumps<-dd$jumps+1
        val$jumpstime<-val$time[val$jumps]
        #val$nevent<-length(val$S0)
        return(val)
      }
      with(val, structure(nevent,gradient=-gradient, hessian=-hessian))
    }
  }
  
  opt <- NULL
  if (p>0) {
    opt<-lava::NR(beta, obj,...)
    opt$estimate <- opt$par
    cc <- opt$estimate; names(cc) <-colnames(X)
    if (!stderr) return(cc)
    val <-c(list(coef=cc), obj(opt$estimate, all=TRUE))
  } else {
    val <- obj(0,all=TRUE)
    val[c("gradient","hessian","U")]<-NULL
  }
  
  res <- c(val,
           list(strata=strata,
                entry=entry,
                exit=exit,
                status=status,
                p=p,
                X=X,
                weight=weight, opt=opt,
                nstrata=nstrata,
                strata.name=strata.name))
  
  class(res) <- "matchpropexc"
  res
}
####}}} matchpropexc0

####{{{ matchpropexc
##' Excess risk paired survival model
##' @param formula formula with 'Surv' outcome (see \code{coxph}); use strata() for strata-variables and cluster() to specify the variable cluster
##' @param data data frame - already set-up. (see\code{compdata} for modre details)
##' @param cluster vector case indicator - by default=cluster
##' @param idControl vector control indicator - by default=unexp.subj
##' @param weight vector that counts how many time the exposed had the event before unexposed invdividuals - by default weight
##' @examples 
##' dd<-data.sim(nca=5000, ncont=5)
##' setdd<-compdata(Surv(time, status)~x+z+cc, cluster=id, idControl=j, data=dd)
##' names(setdd) #it is strongly recommended to check the names of your variables before estimating the model
##' exc.model<-matchpropexc(Surv(exit,status)~strata(z)+factor(x), data=sdd, weight=weight, idControl=unexp.subj, cluster=cluster)
##' summary(exc.model)
##' exc.model1<-matchpropexc(Surv(exit,status)~1, data=sdd,weight=weight, idCluster=unexp.subj, cluster=cluster)
##' summary(exc.model1)
##' @return no output. use \code{summary(model)} to view the coefficient estimates.
##' @author Cristina Boschini
##' @export
matchpropexc <- function(formula, data, cluster, idControl, weight,...){
  if (missing(cluster)) stop("cluster vector needed - use cluster in compdata results")
  if (missing(idControl)) stop("idControl vector needed - use unexp.subj in compdata results")
  if (missing(weight)) stop("cluster weight needed - use weight in compdata results")
  
  cl <- match.call()
  m <- match.call(expand.dots=TRUE)[1:6]
  special <- c("strata")
  Terms <- terms(formula,special,data=data, idControl=idControl, cluster=cluster, weight=weight)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m,"response")
  if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
  if (ncol(Y)==2) {
    exit <- Y[,1]
    entry <- rep(0,nrow(Y))
    status <- Y[,2]
    Truncation <- FALSE
  } else {
    entry <- Y[,1]
    exit <- Y[,2]
    status <- Y[,3]
    Truncation <- TRUE
  }
  strata <- NULL
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)){
    ts <- survival::untangle.specials(Terms, "strata")
    Terms <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name<-ts$vars
  } else strata.name <- NULL
  # cluster <- NULL
  # if (!is.null(attributes(Terms)$specials$cluster)){
  #   ts <- survival::untangle.specials(Terms, "cluster")
  #   Terms <- Terms[-ts$terms]
  #   cluster <- m[[ts$vars]]
  # }
  
  cluster<-model.extract(m,"cluster")
  names(cluster)<- NULL
  idControl<-model.extract(m,"idControl")
  names(idControl)<- NULL
  weight<-model.extract(m,"weight")
  names(weight)<- NULL
  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  p<-ncol(X)
  if(!is.null(colnames(X))) namesX<-colnames(X)
  else if(p>0) namesX <- paste("var",seq(1,p),sep="")
  
  res <- c(matchpropexc0(X,entry=entry,exit=exit,
                         status=status,weight=weight,
                         strata=strata,
                         strata.name=strata.name,...),
           list(call=cl, model.frame=m))
  class(res) <- "matchpropexc"
  
  res
  
}
###}}} matchpropexc

###{{{ sandEst
sandEst<- function(x,...){
  invhess <- solve(x$hessian)
  if (!is.list(x$weight)){ 
    indicator <- ifelse(x$weight>0,1,2)
  } else indicator <- ifelse(unlist(x$weight)>0,1,2)
  ii <- mets::cluster.index(indicator)
  UU <- matrix(nrow=ii$uniqueclust, ncol=ncol(invhess))
  UUpos <- x$U[ii$idclustmat[1,seq(ii$cluster.size[1])]+1,,drop=FALSE]
  UUneg <- x$U[ii$idclustmat[2,seq(ii$cluster.size[2])]+1,,drop=FALSE]
  res <- invhess%*%(crossprod(UUpos)+crossprod(UUneg))%*%invhess
  structure(res, invhess=invhess)
}
###}}} sandEst

###{{{ vcov
##' variance and covariance matrix for the coefficient estimates 
##' @param object model estimated with matchpropexc
##' @author Cristina Boschini
##' @return variance and covariance matrix of the coefficient estimates
##' @export
vcov.matchpropexc <- function(object,...){
  res <-  sandEst(object)
  attributes(res)$invhess <- attributes(sandEst)$invhess
  colnames(res) <- rownames(res) <- names(coef(object))
  res
}
###}}} vcov

###{{{ coef
##' coefficient estimates
##' @param object model estimated with matchpropexc
##' @author Cristina Boschini
##' @export
coef.matchpropexc <- function(object,...) {
  object$coef
}
###}}} coef


###{{{ summary
##' model summary 
##' @param object model estimated with matchpropexc
##' @author Cristina Boschini
##' @return coefficient estimates, standard error and significance level
##' @export
summary.matchpropexc <- function(object,...){
  cc <- NULL
  if (object$p>0) {
    V <- vcov(object)
    cc <- cbind(coef(object),diag(V)^0.5)
    cc <- cbind(cc,2*(pnorm(abs(cc[,1]/cc[,2]), lower.tail=FALSE)))
    colnames(cc) <- c("Estimate","S.E.","P-value")
    rownames(cc) <- names(coef(object))
  }
  Strata<-levels(object$strata)
  if (!is.null(Strata)) {
    n<-unlist(lapply(object$time, length))
  } else {
    n<-length(object$time)
  }
  res <- list(coef=cc,n=n,nevent=object$nevent,
              strata=Strata)
  class(res) <- "summary.matchpropexc"
  res
}

###}}} summary

###{{{ print.summary

##' @export
print.summary.matchpropexc <- function(x,max.strata=5,...){
  cat("\n")
  nn <- cbind(x$n,x$nevent)
  rownames(nn) <- levels(x$strata); colnames(nn) <-c("n","events")
  if (is.null(rownames(nn))) rownames(nn)<-rep("",NROW(nn))
  if (length(x$strata)> max.strata) {
    nn <-rbind(c(colSums(nn),length(x$strata)));
    colnames(nn)<-c("n","events","stratas")
    rownames(nn)<-""
  }
  print(nn,quote=FALSE)
  if (!is.null(x$coef)) {
    cat("\n")
    printCoefmat(x$coef,...)
  }
}

###}}} print.summary

# cumsumstratamc<-function(x,strata,nstrata)
# { # {{{
#   res<-.Call("cumsumstrataR",x,strata,nstrata)$res
#   return(res)
# }}}}


### {{{ vcovCH.mc
vcovCH.mc<-function(p, weight, nevent, X, E, S0, sigmaH, hessian){
  if (p>0){
    l<-2*p
    H<--1*E*c(weight)/c(S0)
    covT <- (X-E)*c(weight^2)/c(S0)
    I <- solve(hessian)
    cumH<-apply(H,2,cumsum)
    
    term1 <- apply(cumH,1,function(x) t(x)%*%sigmaH%*%x)
    term2term3 <- cumsum(weight^2/S0^2)
    term4term5 <- apply(cbind(cumH,apply(covT,2,cumsum)),1,function(x) {
      x1<-x[1:p]
      x2<-x[(p+1):l]
      2*t(x1)%*%I%*%x2
    })
  } else {
    term1<-0
    term2term3<-cumsum(weight^2/S0^2)
    term4term5<-0
  }
  
  vcovchaz <- term1+term2term3+term4term5
  return(vcovchaz)
}

### }}} vcovCH.mc

### {{{ exccumhaz

cumhazmc<-function(time, weight, S0, p, nevent, X, E, sigmaH=NULL, hessian, SEcumhaz=TRUE) {
  if (SEcumhaz) {
    chaz <- cbind(time, cumsum(weight/S0),(vcovCH.mc(p,weight, nevent, X, E, S0, sigmaH,hessian))^0.5)
    colnames(chaz)<-c("time","chaz","se.chaz")
    
  }
  else {
    chaz <- cbind(time, cumsum(weight/S0))
    colnames(chaz)<-c("time","chaz")
  }
  return(chaz)
}

##' Cumulative baseline excess hazard estimate
##' @param object model estimated with matchpropexc
##' @param strata default strata if the model has them. 
##' @param time specify at which time to compute the estimates
##' @param SEcumhaz by default set to TRUE. FALSE if you don't want to compute it.
##' @author Cristina Boschini
##' @return cumulative baseline excess hazard estimates. If the model has strata, the returned object will be a list. Estimated are computed at the defined time.
##' @export
exccumhaz<-function(object, strata=object$strata, time=NULL,
                        SEcumhaz=TRUE){
  if (object$p>0) sigmaH <- vcov(object)
  if (is.null(strata)) {
    chaztab<-cumhazmc(object$jumpstime, object$weight, object$S0, 
                      object$p, object$nevent, object$xjumps, object$E,
                      sigmaH, object$hessian, SEcumhaz)
    if (!is.null(time)) {
      if (SEcumhaz) {
        chaztab<-cbind(timereg::Cpred(chaztab[,1:2], time),
                       timereg::Cpred(chaztab[,c(1,3)], time)[,-1])
        colnames(chaztab)<-c("time","chaz","se.chaz")
      } else {
        chaztab<-timereg::Cpred(chaztab, time)
        colnames(chaztab)<-c("time","chaz")
      }
    }
  } else {  
    lev<-levels(strata)
    chaztab<-c()
    for (i in seq(length(lev))){
      chaztab<-c(chaztab, list(cumhazmc(object$jumpstime[[i]], object$weight[[i]], object$S0[[i]], 
                        object$p, object$nevent[[i]], object$xjumps[[i]], object$E[[i]],
                        sigmaH, object$hessian, SEcumhaz)))
    }
    names(chaztab)<-lev
    if(!is.null(time)){
      if (SEcumhaz) {
        chaztab<-lapply(chaztab, function(x) {
          tab<-cbind(timereg::Cpred(x[,1:2], time),
                     timereg::Cpred(x[,c(1,3)], time)[,-1])
          colnames(tab)<-c("time","chaz","se.chaz")
          return(tab)
      }
      )
      } else {
        chaztab<-lapply(chaztab, function(x) {
          tab<-timereg::Cpred(x[,1:2], time)
          colnames(tab)<-c("time","chaz")
          return(tab)
        }
        )
      }
    }
  }
    return(chaztab)
}

### }}} cumhaz


###{{{ predict 

predictmc<- function(chaztab, beta,X=NULL,relsurv=FALSE,...){
  if (!ncol(chaztab)==3) warning("no SE for cumhaz")
  chaz<-chaztab

  if( !is.null(X)) {
    H<-exp(X%*%beta)
    if (nrow(chaz)==length(H)) {
      chaz <- cbind(chaz[,1],chaz[,2]*H)
    } else {
      chaz2 <- c()
      X <- rbind(X)
      for (i in seq(nrow(X)))
        chaz2 <- rbind(chaz2,
                       cbind(chaz[,1],chaz[,2]*H[i],
                             rep(1,nrow(chaz))%x%X[i,,drop=FALSE]))
      chaz <- chaz2;
      nn <- c("time","chaz",names(beta))
      colnames(chaz) <- nn
    }
  }
  if (relsurv) {
    chaz[,2]<-exp(-chaz[,2])
    colnames(chaz)[2]<-"relsurv"
  }
  return(chaz)
}

##' Predict cumulative excess hazard or relative survival
##' @param object model estimated with matchpropexc
##' @param relsurv set to FALSE. If TRUE relative survival is returned instead of the cumulative baseline hazard
##' @param time specify at which time to compute the estimates
##' @param X define specific values for the excess covariates. By default predictworks on the original data.
##' @param strata default strata if the model has them.
##' @author Cristina Boschini
##' @export
predict.matchpropexc <- function(object, 
                                 relsurv=FALSE,
                                 time=object$exit,
                                 X=object$X,
                                 strata=object$strata,...){
  if (object$p==0) X<-NULL
  if(!is.null(object$strata)){
    lev <-levels(object$strata)
    if (!is.null(object$strata) &&
        !(is.list(time) & !is.data.frame(time)) &&
        !(is.list(X) & !is.data.frame(X))) {
      X0 <- X
      time0 <- time
      X <- time <- c()
      if (all(time0 %in% object$exit)){
        for (i in seq(length(lev))) {
          idx <- which(strata==lev[i])
          X <- c(X,list(X0[idx,,drop=FALSE]))
          time <- c(time,list(time0[idx]))
        }
      } else {
        for (i in seq(length(lev))) {
          idx <- which(strata==lev[i])
          X <- c(X,list(X0[idx,,drop=FALSE]))
          time<-c(time,list(time0))
        }
      }
    }
    chaz<-c()
    for (i in seq(length(lev))) {
      chaztab<-exccumhaz(object,time=time[[i]])[[i]]
      chaz<-c(chaz, list(predictmc(chaztab,
                                   coef(object),
                                   X[[i]],relsurv)))
    }
    names(chaz)<-lev
  } else {
    chaztab<-exccumhaz(object,time=time)
    chaz <- predictmc(chaztab,coef(object),X, relsurv)
  }
  return(chaz)
}
###}}} predict

###{{{ excplot

##' Plotting the cumulative baseline excess hazard
##'
##' Plotting the cumulative baselines excess hazard
##' @param x matchpropexc object
##' @param se to include standard errors
##' @param time to plot for specific time variables
##' @param add to add to previous plot 
##' @param ylim to give ylim 
##' @param lty to specify lty of components
##' @param col to specify col of components
##' @param legend to specify col of components
##' @param ylab to specify ylab 
##' @param polygon to get standard error in shaded form
##' @param level of standard errors
##' @param stratas wich strata to plot (number or vector between 0 and nstrata-1)
##' @param relsurv set to FALSE. If TRUE relative survival curves are plotted.
##' @param ... Additional arguments to lower level funtions
##' @author Cristina 
##' @examples 
##' excplot(m, se=TRUE, col=c("green","blue")) #with condifence intervals
##' excplot(m, se=TRUE, time=seq(0,30,1)) #at specific time-points
##' excplot(m, se=TRUE, relsurv=TRUE) #plot the relative survival
##' excplot(m, se=TRUE, polygon=FALSE) #No polygon
##' excplot(m, se=FALSE) #no confidence interval
##' excplot(m, se=TRUE, stratas=1) #plot the second strata
##' @export
excplot  <- function(x, se=FALSE,
                                  time=NULL, add=FALSE,
                                  ylim=NULL,
                                  lty=NULL,col=NULL,legend=TRUE,
                                  ylab="Baseline cumulative excess hazard",
                                  polygon=TRUE,
                                  level=0.95,
                                  stratas=NULL,
                                  relsurv=FALSE,...) {# {{{
  level <- -qnorm((1-level)/2)
  ## all strata
  if (is.null(stratas)) stratas <- 0:(x$nstrata-1) 
  if (!is.null(stratas) & !any(stratas %in% 0:(x$nstrata-1))) stop("The selected strata doesn't exist\n")
  
  ltys <- lty
  cols <- col
  

  ## with strata
  if (length(stratas)>0 & x$nstrata>1) { 
    ms <- match(x$strata.name,names(x$model.frame))
    lstrata <- levels(x$model.frame[,ms])[(stratas+1)]
    stratn <-  substring(x$strata.name,8,nchar(x$strata.name)-1)
    stratnames <- paste(stratn,lstrata,sep=":")
    if (!is.matrix(lty)) {
      if (is.null(lty)) ltys <- 1:length(stratas) else if (length(lty)!=length(stratas)) ltys <- rep(lty[1],length(stratas))
    } else ltys <- lty
    if (!is.matrix(col)) {
      if (is.null(col)) cols <- 1:length(stratas) else 
        if (length(col)!=length(stratas)) cols <- rep(col[1],length(stratas))
    } else cols <- col
  } else { 
    stratnames <- "Baseline" 
    if (is.matrix(col))  cols <- col
    if (is.null(col)) cols <- 1  else cols <- col[1]
    if (is.matrix(lty))  ltys <- lty
    if (is.null(lty)) ltys <- 1  else ltys <- lty[1]
  }
  
  if (!is.matrix(ltys))  ltys <- cbind(ltys,ltys,ltys)
  if (!is.matrix(cols))  cols <- cbind(cols,cols,cols)
  
  listcumhaz<-exccumhaz(x, time=time)
  if(!is.null(x$strata)) {
    cumhaz<-do.call("rbind",listcumhaz )
  } else cumhaz<-listcumhaz
  
  if(!relsurv) {
    rr <- range(cumhaz[,2])
  } else rr<-range(exp(-cumhaz[,2]))
  
  if (is.null(ylim)) ylim <- rr
  if (se==TRUE) {
    if (ncol(cumhaz)<3) stop("exccumhaz must be with SE.cumhaz=TRUE\n"); 
    if (!relsurv) {
      rrse <- range(c(cumhaz[,2]+level*cumhaz[,3]))
    } else {
      rrse <- range(c(exp(-(cumhaz[,2]+level*cumhaz[,3]))))
    }
    ylim <- rrse
  }

  i<-1  
  
  if(!is.null(x$strata)) {
    if(length(stratas)>1) {
      listcumhaz1<-listcumhaz[[i+stratas[1]]]
    } else listcumhaz1<-listcumhaz[[i+stratas]]
  } else listcumhaz1<-listcumhaz
  
  if (!relsurv) {
    plotcurve <- listcumhaz1[,1:2]
  } else {
    plotcurve<-cbind(listcumhaz1[,1], exp(-listcumhaz1[,2]))
  }
  
  cumhazse<-listcumhaz1
  
  if (add) {
    lines(plotcurve,type="s",lty=ltys[i,1],col=cols[i,1],...)
  } else {
    plot(plotcurve,type="s",lty=ltys[i,1],col=cols[i,1],ylim=ylim,ylab=ylab,xlab="Time",...)
  }
  if (se==TRUE) {
    if (!relsurv) {
      secumhazard <- cumhazse[,c(1,3)]
      ul <-cbind(plotcurve[,1],plotcurve[,2]+level*secumhazard[,2])
      nl <-cbind(plotcurve[,1],plotcurve[,2]-level*secumhazard[,2])
    } else {
      secumhazard <- cumhazse[,c(1,3)]
      cumhazard <- cumhazse[,2]
      ul <-cbind(plotcurve[,1],exp(-(cumhazard-level*secumhazard[,2])))
      nl <-cbind(plotcurve[,1],exp(-(cumhazard+level*secumhazard[,2])))
    }
    
    if (!polygon) {
      lines(nl,type="s",lty=ltys[i,2],col=cols[i,2])
      lines(ul,type="s",lty=ltys[i,3],col=cols[i,3])
    } else {
      tt <- c(nl[,1],rev(ul[,1]))
      yy <- c(nl[,2],rev(ul[,2]))
      col.alpha<-0.1
      col.ci<-cols[i]
      col.trans <- sapply(col.ci, FUN=function(x) 
        do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))
      polygon(tt,yy,lty=ltys[i,2],col=col.trans, border = cols[i,1])
    }
  }
  
  if (length(stratas)>1)  {
    for (i in 2:length(stratas)) {
      if (!relsurv) {
        plotcurve <- listcumhaz[[i]][,1:2]
      } else {
        plotcurve<-cbind(listcumhaz[[i]][,1], exp(-listcumhaz[[i]][,2]))
      }
      cumhazse<-listcumhaz[[i]]
      lines(plotcurve,type="s",lty=ltys[i,1],col=cols[i,1])   
      if (se==TRUE) {
        if (!relsurv) {
          secumhazard <- cumhazse[,c(1,3)]
          ul <-cbind(plotcurve[,1],plotcurve[,2]+level*secumhazard[,2])
          nl <-cbind(plotcurve[,1],plotcurve[,2]-level*secumhazard[,2])
        } else {
          secumhazard <- cumhazse[,c(1,3)]
          cumhazard<-cumhazse[,2]
          ul <-cbind(plotcurve[,1],exp(-(cumhazard-level*secumhazard[,2])))
          nl <-cbind(plotcurve[,1],exp(-(cumhazard+level*secumhazard[,2])))
        }
        
        if (!polygon) {
          lines(nl,type="s",lty=ltys[i,2],col=cols[i,2])
          lines(ul,type="s",lty=ltys[i,3],col=cols[i,3])
        } else {
          tt <- c(nl[,1],rev(ul[,1]))
          yy <- c(nl[,2],rev(ul[,2]))
          col.alpha<-0.1
          col.ci<-cols[i]
          col.trans <- sapply(col.ci, FUN=function(x) 
            do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))
          polygon(tt,yy,lty=ltys[1,2],col=col.trans, border = cols[i,1])
        }
      }
    }
  }
  
  if (legend)
    if(!relsurv) {
      legend("topleft",legend=stratnames,col=cols[,1],lty=ltys[,1], bty="n")
    } else legend("bottomleft",legend=stratnames,col=cols[,1],lty=ltys[,1], bty="n")
  
  
}# }}} 

##' @export
lines.matchpropexc <- function(x,...,add=TRUE) excplot(x,...,add=add)

###}}} excplot

###{{{ plot

##' @export
plot.matchpropexc  <- function(x,relsurv=TRUE,X=NULL,time=NULL,add=FALSE,...) {
  browser()
  if (!is.null(X) && nrow(X)>1) {
    P <- lapply(split(X,seq(nrow(X))),function(xx) predict(x,X=xx,time=time,surv=surv))
  } else {
    P <- predict(x,X=X,time=time,surv=surv)
  }
  if (!is.list(P)) {
    if (add) {
      lines(P,type="s",...)
    } else {
      plot(P,type="s",...)
    }
    return(invisible(P))
  }
  
  if (add) {
    lines(P[[1]][,1:2],type="s",lty=1,col=1,...)
  } else {
    plot((P[[1]])[,1:2],type="s",lty=1,col=1,...)
  }
  for (i in seq_len(length(P)-1)+1) {
    lines(P[[i]][,1:2],type="s",lty=i,col=i,...)   
  }
  return(invisible(P))
}

##' @export
lines.phreg <- function(x,...,add=TRUE) plot(x,...,add=add)

###}}} plot


