### {{{ compdata
compdata<-function(entry,exit,status,cluster,idControl,X,strata,Truncation){
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
    stat <- ifelse(indicator==-1,1,ifelse(weight==1,0,ifelse(weight==0,0,1)))
    
    d3 <- data.frame(start,end,indicator,stat,group,subj,weight)
    
    if (ncol(X)>0) {
      Xcases <- as.matrix(X[wp$pairs[,2],], ncol=ncol(X))
      colnames(Xcases)<-colnames(X)
      d3 <- data.frame(d3,Xcases, check.names=FALSE)
    }
    
    if (!is.null(strata)) {
      strata <- strata[wp$pairs[,2]]
      d3 <- data.frame(d3, strata, check.names = FALSE)
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
    stat <- ifelse(indicator==-1,1,ifelse(weight==1,0,ifelse(weight==0,0,1)))
    
    d3 <- data.frame(end,indicator,stat,group,subj,weight)
    
    if (ncol(X)>0) {
      Xcases <- as.matrix(X[wp$pairs[,2],], ncol=ncol(X))
      colnames(Xcases)<-colnames(X)
      d3 <- data.frame(d3,Xcases, check.names=FALSE)
    }
    
    if (!is.null(strata)) {
      strata <- strata[wp$pairs[,2]]
      d3 <- data.frame(d3, strata, check.names = FALSE)
    } 
  }
  return(d3)
  
}

###}}} compdata

###{{{ matchpropexc0

matchpropexc0 <- function(X,entry, exit, status, weight,strata=NULL, beta,stderr=TRUE,cumhaz=TRUE,...){
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
                weight=weight, opt=opt))
  
  class(res) <- "matchpropexc"
  res
}
####}}} matchpropexc0

####{{{ matchpropexc
##' Excess risk paired survival model
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param idControl vector control indicator (the number of controls is expected to be the same for every case)
##' @author Cristina Boschini
##' @export
matchpropexc <- function(formula,data,idControl,...){
  #  idCase <- eval(substitute(idCase),data)
  #  idControl <- eval(substitute(idControl),data)
  #  if (is.null(idCase) | is.null(idControl)) stop("idCase and idControl needed")
  cl <- match.call()
  m <- match.call(expand.dots=TRUE)[1:4]
  special <- c("strata","cluster")
  Terms <- terms(formula,special,data=data, idControl=idControl)
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
  }
  cluster <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)){
    ts <- survival::untangle.specials(Terms, "cluster")
    Terms <- Terms[-ts$terms]
    cluster <- m[[ts$vars]]
  }
  idControl<-model.extract(m,"idControl")
  names(idControl)<- NULL
  X <- model.matrix(Terms, m)
  if (!is.null(intpos <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  p<-ncol(X)
  if(!is.null(colnames(X))) namesX<-colnames(X)
  else if(p>0) namesX <- paste("var",seq(1,p),sep="")
  
  if (Truncation) {
    setupdata <- compdata(entry,exit,status,cluster,idControl,X,strata,Truncation)
  } else setupdata <- compdata(entry=NULL,exit,status,cluster,idControl,X,strata,Truncation)
  
  exit <-setupdata$end    
  if (Truncation) {
    entry <- setupdata$start
  } else entry <- rep(0,nrow(Y))
  
  if (p>0) {
    X <- as.matrix(setupdata[,namesX], ncol=p)
    colnames(X)<-namesX
  }
  
  status <- setupdata$stat
  weight <- setupdata$weight
  strata <- setupdata$strata
  res <- c(matchpropexc0(X,entry, exit,status,weight,strata,...),list(call=cl, model.frame=m))
  class(res) <- "matchpropexc"
  
  res
  
}
###}}} matchpropexc

###{{{ vcov
##' @export
vcov.matchpropexc <- function(object,...){
  res <-  sandEst(object)
  attributes(res)$invhess <- attributes(sandEst)$invhess
  colnames(res) <- rownames(res) <- names(coef(object))
  res
}
###}}} vcov

###{{{ coef
##' @export
coef.matchpropexc <- function(object,...) {
  object$coef
}
###}}} coef

###{{{ sandEst
##' @export
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

###{{{ summary

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
vcovCH.mc<-function(p, weight, nevent, X, E, S0, sigmaH=NULL, hessian){
  if (p>0){
    l<-2*p
    
    matW <- diag(as.vector(weight), nrow=nevent)
    H <- apply((matW%*%E), 2, function(x) (x/S0)*(-1))
    covT <- apply((X-E),2, function(x) (x/S0)*weight^2)
    I <- solve(hessian)
    
    term1 <- apply(apply(H,2,cumsum),1,function(x) t(x)%*%sigmaH%*%x)
    term2term3 <- cumsum(weight^2/S0^2)
    term4term5 <- apply(cbind(apply(H,2,cumsum),apply(covT,2,cumsum)),1,function(x) {
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

### {{{ cumhaz.matchf
##' @export
cumhaz.matchf<-function(object, strata=object$strata){
  browser() 
  if (object$p>0) sigmaH <- vcov(object)
  if (is.null(strata)) {
  ### computes Breslow estimator and SE
  jumpstime<-object$jumpstime
  chaz <- cbind(jumpstime, cumsum(object$weight/object$S0))
  se.chaz<-cbind(jumpstime,(vcovCH.mc(object$p,object$weight, object$nevent, object$xjumps, object$E, object$S0, sigmaH,
                                     object$hessian))^0.5)
  colnames(chaz)<-c("time","chaz")
  colnames(se.chaz)<-c("time","se.chaz")
  } else {  
    lev<-levels(strata)
    chaz<-c()
    se.chaz<-c()
    for (i in seq(length(lev)))
      chaz<-c(chaz, list(cbind(object$jumpstime[[i]],
                              cumsum(object$weight[[i]]/object$S0[[i]]))))
      se.chaz<-c(se.chaz, list(vcovCH.mc(object$p,
                                         object$weight[[i]],
                                         object$nevent[[i]],
                                         object$xjumps[[i]],
                                         object$E[[i]],
                                         object$S0[[i]],
                                         object$hessian,
                                         sigmaH)))
      names(chaz)<-names(se.chaz)<-lev
  }
}

### }}} cumhaz.matchf


###{{{ predict

predictmc<- function(x,jumpstime, S0, weight, beta, time=NULL,X=NULL,relsurv=FALSE,...){
  if( !is.null(X)) {
    H<-exp(X%*%beta)
    if (nrow(chaz)==length(H)) {
      chaz[,2] <- chaz[,2]*H
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

##' @export
predict.matchpropexc <- function(object, data,
                                 time=object$exit,
                                 X=object$X,
                                 strata=object$strata,
                                 relsurv=FALSE,...){
  if (object$p==0) X<-NULL
  if (!is.null(object$strata) &&
      !all(time %in% object$exit) &&
      !(is.list(X) & !is.data.frame(X))) {
    lev <-levels(object$strata)
    time0<-time
    time<-rep(list(time0), length(lev))
    X0<-X
    X<-rep(list(X), length(X))
  }
  if(!is.null(object$strata)) {
    lev <-levels(object$strata)
    if (!is.null(object$strata) &&
        !(is.list(time) & !is.data.frame(time)) &&
        !(is.list(X) & !is.data.frame(X))) {
      time0<-time
      X0<-X
      X<-time<-c()
      for (i in seq(length(lev))) {
        idx <- which(strata==lev[i])
        X <-c(X, list(X0[idx,,drop=FALSE]))
        time <-c(time, list(time0[idx]))
      }
    }
    chaz<-c()
    for (i in seq(length(lev)))
      chaz<-c(chaz, list(predictmc(object$jumpstime[[i]],
                                   object$S0[[i]],
                                   object$weight[[i]],
                                   coef(object),
                                   time[[i]],
                                   X[[i]], relsurv)))
    names(chaz)<-lev
  } else {
    chaz <- predictmc(object$jumpstime, object$S0, object$weight, coef(object),time,
                      X, relsurv)
  }
  return(chaz)
}
###}}} predict

###{{{ plot

##' @export
plot.matchpropexc <-function(x, relsurv=FALSE, X=NULL, time=NULL, add=FALSE, conf.int=FALSE,...) {
  require(ggplot2)
  if (!is.null(X) && nrow(X)>1) {
    P <- lapply(unique(split(X,seq(nrow(X)))),function(xx) predict(x,X=xx,time=time,relsurv=surv))
  } else {
    P <- predict(x,X=X,time=time,relsurv=relsurv)
  }
  
  
}
###}}}

