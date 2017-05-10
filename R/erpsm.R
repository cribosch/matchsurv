### {{{ compdata
compdata<-function(entry,exit,status,cluster,idControl,X,strata,Truncation){
    if (Truncation){
       
        ii <- mets::cluster.index(cluster)
        jj <- mets::cluster.index(idControl)
        
        k <- ii$maxclust-1
        
        matentry <- t(matrix(entry[jj$idclustmat+1], ncol=jj$maxclust))
        matexit <- t(matrix(exit[jj$idclustmat+1], ncol=jj$maxclust))
        matstatus <- t(matrix(status[jj$idclustmat+1], ncol=jj$maxclust))
	if (ncol(X)>0) {
            if(!is.null(colnames(X))) {namesX<-colnames(X)
            } else namesX <- paste("var", seq(1,ncol(X)),sep="")
        }
        
        pd2<-NULL
        for (j in 1:k){
            EntryMax <- pmax(matentry[,1],matentry[,j+1], na.rm=TRUE)
            EntryMax <- ifelse(is.na(matentry[,j+1]),NA, EntryMax)
            ExitMin <- pmin(matexit[,1],matexit[,j+1], na.rm=TRUE)
            ExitMin <- ifelse(is.na(matexit[,j+1]),NA, ExitMin)
            
            Smin<-ifelse(is.na(ExitMin),NA,
                  ifelse(ExitMin<EntryMax,NA,
                  ifelse(ExitMin==matexit[,1] & EntryMax<matexit[,1],
                         matstatus[,1],
                  ifelse(matstatus[,j+1]==1 & EntryMax<=ExitMin,-1,0))))
            
            pd2<-cbind(pd2,cbind(EntryMax, ExitMin,Smin))
        }
        pd2<-data.frame(pd2)    
        pd2$w<-rowSums(pd2[,seq(3,k*3,3)]==1,na.rm=TRUE)
        pd2$id<-1:jj$maxclust
        if (is.vector(X)) X<-matrix(X, ncol=1)

        if (!is.null(strata)) {
            pd2<-data.frame(pd2,strata[!duplicated(idCase)])
        } else {
            pd2<-data.frame(pd2)
        }

        if (ncol(X)>0) {
            Xcase <- X[!duplicated(idCase),]
            pd2 <- data.frame(pd2,Xcase)
        }
    
        d2<-data.table::data.table(pd2)
        
        if (!is.null(strata)) {
            colnamesd2<-c(paste0(c("entry","exit","indicator"),
                                   rep(seq_len(k),each=3)),
                            "w","id", "strata")
        } else {
            colnamesd2<-c(paste0(c("entry","exit","indicator"),
                                   rep(seq_len(k),each=3)),
                            "w","id")
        }

        if (ncol(X)>0) {
            colnamesd2<-c(colnamesd2, namesX)
        }

        colnames(d2) <- colnamesd2
        
        d3<-data.table::melt(d2,measure=patterns("^entry","^exit","^indicator"),
                             value.name=c("entry","exit","indicator"),
                             variable.factor=FALSE, variable.name="j")
        class(d3)<-"data.frame"
        d3$status<-1
        d3$w<-ifelse(d3$indicator==0,1,d3$w)
        d3$status<-ifelse(d3$indicator==0,0,d3$status)
        d3$w<-ifelse(d3$indicator==-1,-1,d3$w)
        d3$dummy<-paste0(d3$id,".",d3$w)
        d3$w <-ifelse(d3$w>1 & duplicated(d3$dummy),0,d3$w)
        d3$status<-ifelse(d3$w==0,0,d3$status)
        d3<-d3[d3$exit>d3$entry,]
    } else {
        
        ii <- mets::cluster.index(idCase)
        jj <- mets::cluster.index(idControl)
        
        k <- ii$maxclust-1
        
        matexit <- t(matrix(exit[jj$idclustmat+1], ncol=jj$maxclust))
        matstatus <- t(matrix(status[jj$idclustmat+1], ncol=jj$maxclust))
        
        if (ncol(X)>0) {
            if(!is.null(colnames(X))) {namesX<-colnames(X)
            } else namesX <- paste("var", seq(1,ncol(X)),sep="")
        }
        

        pd2<-NULL
        for (j in 1:k){
            ExitMin <- pmin(matexit[,1],matexit[,j+1], na.rm=TRUE)
            ExitMin <- ifelse(is.na(matexit[,j+1]),NA, ExitMin)
            
            Smin<-ifelse(is.na(ExitMin),NA,
                  ifelse(ExitMin==matexit[,1], matstatus[,1],
                  ifelse(matstatus[,j+1]==1,-1,0)))
            
            pd2<-cbind(pd2,cbind(ExitMin,Smin))
        }
        
        pd2<-data.frame(pd2)    
        pd2$w<-rowSums(pd2[,seq(2,k*2,2)]==1,na.rm=TRUE)
        pd2$id<-1:jj$maxclust
        if (is.vector(X)) X<-matrix(X, ncol=1)

        if (!is.null(strata)) {
            pd2<-data.frame(pd2, strata[!duplicated(idCase)])
        } else {
            pd2<-data.frame(pd2)
        }
        
        if (ncol(X)>0) {
            Xcase <- X[!duplicated(idCase),]
            pd2 <- data.frame(pd2, Xcase)
        }
        
        d2<-data.table::data.table(pd2)
        
        if (!is.null(strata)) {
            colnamesd2<-c(paste0(c("exit","indicator"),
                               rep(seq_len(k), each=2)),
                            "w","id", "strata")
        } else {
            colnamesd2<-c(paste0(c("exit","indicator"),
                               rep(seq_len(k), each=2)),
                            "w","id")
        }
        
        if (ncol(X)>0) {
            colnamesd2<-c(colnamesd2, namesX)
        }
        colnames(d2) <- colnamesd2
        
        d3<-data.table::melt(d2,measure=patterns("^exit","^indicator"),
                             value.name=c("exit","indicator"),
                             variable.factor=FALSE, variable.name="j")
        
        
        class(d3)<-"data.frame"
        d3$status<-1
        d3$w<-ifelse(d3$indicator==0,1,d3$w)
        d3$status<-ifelse(d3$indicator==0,0,d3$status)
        d3$w<-ifelse(d3$indicator==-1,-1,d3$w)
        d3$dummy<-paste0(d3$id,".",d3$w)
        d3$w <-ifelse(d3$w>1 & duplicated(d3$dummy),0,d3$w)
        d3$status<-ifelse(d3$w==0,0,d3$status)
        
    }
    d3$dummy<-NULL
    return(d3)
    
}

###}}} compdata


###{{{ ersp0

erpsd0 <- function(X,entry, exit, status, weight,strata=NULL, beta,stderr=TRUE,...){
    if(is.vector(X)) X <- matrix(X, ncol=1)
    p <-ncol(X)
    if (missing(beta)) beta <-rep(0,p)
    if (p==0) X<-cbind(rep(0,length(exit)))
    if (!is.null(strata)){
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
                jumptimes<-lapply(dd, function(x) x$time[x$ord+1][x$jumps+1])
                weight<-lapply(val, function(x)  x$weight)
                #S0 <- lapply(val, function(x) x$S0)
                #nevent<-unlist(lapply(S0,length))
                return(list(gradient=gradient, hessian=hessian,
                            U=U, S0=S0, nevent=nevent,
                            ord=ord, time=time, jumps=jumps, 
                            jumptimes=jumptimes, weight=weight))
            }
            structure(nevent,gradient=-gradient, hessian=-hessian)
            }
        } else {

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
            
            class(res) <- "erpsd"
            res
}
####}}} erpsd0

####{{{ erpsd
##' Excess risk paired survival model
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param idCase vector case indicator
##' @param idControl vector control indicator (the number of controls is expected to be the same for every case)
##' @author Cristina Boschini
##' @export
erpsd <- function(formula,data,idControl,...){
	browser()
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
    X <- model.matrix(Terms, m)
    if (!is.null(intpos <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
    p<-ncol(X)
    if(!is.null(colnames(X))) namesX<-colnames(X)
    else namesX <- paste("var",seq(1,ncol(X)),sep="")

    if (Truncation) {
        setupdata <- compdata(entry,exit,status,cluster,idControl,X,strata,Truncation)
    } else setupdata <- compdata(entry=NULL,exit,status,cluster,idControl,X,strata,Truncation)
    
    exit <-setupdata$exit    
    if (Truncation) {
        entry <- setupdata$entry
    } else entry <- rep(0,nrow(Y))

    if (p>0) {
        X <- as.matrix(setupdata[,namesX], ncol=p)
        colnames(X)<-namesX
    }
    
    status <- setupdata$status
    weight <- setupdata$w
    strata <- setupdata$strata
    res <- c(erpsd0(X,entry, exit,status,weight,strata,...),list(call=cl, model.frame=m))
    class(res) <- "erpsd"
    
    res
    
}
###}}} erpsd

###{{{ vcov
##' @export
vcov.erpsd <- function(object,...){
    res <-  sandEst(object)
    attributes(res)$invhess <- attributes(sandEst)$invhess
    colnames(res) <- rownames(res) <- names(coef(object))
    res
}
###}}} vcov

###{{{ coef
##' @export
coef.erpsd <- function(object,...) {
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
summary.erpsd <- function(object,...){
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
    class(res) <- "summary.erpsd"
    res
}
    
###}}} summary

###{{{ print.summary

##' @export
print.summary.erpsd <- function(x,max.strata=5,...){
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

###{{{ predict

predictErpsd <- function(jumptime, S0, weight, beta, time=NULL,...){
    ## Brewslow estimator
    chaz <- cbind(jumptime, cumsum(weight/S0))
    if (!is.null(time)){
        chaz <- timereg::Cpred(chaz, time)
    }
    colnames(chaz)<-c("time","chaz")

    return(chaz)
}

##' @export
predict.erpsd <- function(object, data, time=object$exit,strata=object$strata,...){
	if (!is.null(strata) && 
	    !any(time %in% object$exit)) {
		lev <-levels(object$strata)
		time0<-time
		time<-rep(list(time0), length(lev))
			  }
	if(!is.null(object$strata)) {
		lev <-levels(object$strata)
		if (!is.null(object$strata) && 
		    !(is.list(time) & !is.data.frame(time))) {
			time0<-time
			time<-c()
			for (i in seq(length(lev))) {
				idx <- which(strata==lev[i])
				time <-c(time, list(time0[idx]))
				}
			} 
		chaz<-c()
		for (i in seq(length(lev)))
			chaz<-c(chaz, list(predictErpsd(object$jumptime[[i]],
							object$S0[[i]],
							object$weight[[i]],
							coef(object),
							time[[i]])))
		names(chaz)<-lev
		} else {
    chaz <- predictErpsd(object$jumpstime, object$S0, object$weight, coef(object),time)
		}
    return(chaz)
}
###}}} predict

###{{{ vcovCH

vcovCHErpsd <- function(jumptime, S0, weight, beta,X,E,hessian,sigmaH,time=NULL,...){
    ## SE for Breslow estimator
    nevent <- length(jumptime)
    p<-ncol(X)
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
    vcovchaz <- cbind(jumptime, term1+term2term3+term4term5)

    if (!is.null(time)){
        vcovchaz <- timereg::Cpred(vcovchaz, time)
    }
    colnames(vcovchaz)<-c("time","vcov.chaz")

    return(vcovchaz)
}

##' @export
vcovCH.erpsd <- function(object, time=object$exit){
    sigmaH <- vcov(object)
    vcovchaz <- vcovCHErpsd(object$jumpstime, object$S0, object$weight, coef(object),
                            object$xjumps,object$E,object$hessian,sigmaH, time)
    return(vcovchaz)
}
###}}} vcovCH
