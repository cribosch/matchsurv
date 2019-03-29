### This script contains the function to estimate the model. No simulation functions. 

### compcomp function: the aim is to get a dataset that then can be used with prep.glm.comprisk in order to estimate the excess cumulative 
### incidence function in a matched cohort.
### It is similar to compdata but for CIF


#globalVariables
if(getRversion() >= "2.15.1")  utils::globalVariables(c("V1",
                                                        "uexit","ucause",
                                                        "eexit","ecause", 
                                                        "clust.num","weightedo",
                                                        "Rt","subj","h",
                                                        "pexittime", "pcause",
                                                        "weights_","cw_","cw"))

### {{{ compcomp
##' Data structured for glm approach in competing risk setting
##' @param formula formula with 'Event' outcome (see \code{timereg} package); time stands for the start time, while time2 stands for the stop time. cause=1 will be considered as the event of interest
##' @param data data frame
##' @param idControl vector control indicator (idControl==1 indicates exposed individual in cluster i)
##' @param cluster vector cluster indicator (one cluster for each exposed individual)
##' @param strata weights computed according to strata
##' @param time.points vector of time points where the glm will be estimated (10 usually is a sufficient number; the more time points, the slower the glm function)
##' @param cens.formula useful to estimate the weights when censoring is present. no quotes, add something like ~age+year
##' @param cens.code default is 0
##' @param event in which event are you interested in?
##' @import data.table
##' @import stats
##' @import survival
##' @examples
##' dcif<-sim.data.MatchCR(nca=1000, ncont=5)
##' tp<-c(0.5,1,2,5,10,15,25)
##' setdcif1<-compcomp(timereg::Event(time=FALSE,time2=time,cause=cause)~X1+X2,
##'  data=dcif, cluster=i, idControl=j, time.points=tp, cens.formula=NULL, event=1)
##' head(setdcif1, 10)
##' setdcif2<-compcomp(timereg::Event(time=FALSE,time2=time,cause=cause)~X1+X2,
##'  data=dcif, cluster=i, idControl=j, time.points=tp, cens.formula=NULL, event=2)
##' head(setdcif2, 10)
##' @author Cristina Boschini
##' @return A setup dataset, ready for \code{geese}
##' @export
compcomp<-function(formula,data,cluster,idControl, strata=NULL,
                   time.points,cens.formula=NULL, cens.code=0, event=1){
  #browser()
  currentOPTs <- options("na.action")
  options(na.action = "na.pass")
  m <- match.call()[1:5]
  Terms <- terms(formula,data=data,cluster=cluster, idControl=idControl, strata=strata)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  cause <- model.response(m)[,3]
  exit <- model.response(m)[,2]
  entry <- model.response(m)[,1]
  if (any(entry<0)) stop("entry time must be positive")
  idControl<-model.extract(m,"idControl")
  names(idControl)<- NULL
  cluster<-model.extract(m,"cluster")
  names(cluster)<- NULL
  strata<-model.extract(m,"strata")
  names(strata)<-NULL
  
  if (sum(entry)>0) {
    Truncation <- TRUE
  } else {
    Truncation <- FALSE
  }
  if (any(entry==exit)) { #I don't know if it is a valid solution. check it!
      warning("Time to event might be null\n")
      entryna<-entry
      exitna<-exit
      exitna[entry==exit]<-(exit[entry==exit]+runif(1,0,0.01))
      entry<-entryna
      exit<-exitna
    }
  
  cord<-order(cluster,idControl)
  
  if(data.table::is.data.table(data)) {
    X <- data[,attributes(Terms)$term.labels, with=FALSE]
    #X<-cbind(X,cluster)
    X<-X[cord]
  } else {
    X<-data[,attributes(Terms)$term.labels, drop=FALSE]
    #X<-cbind(X,cluster)
    X<-X[cord,]
  }
  if (!is.null(strata)) X<-cbind(X,strata)
  if (ncol(X)!=0) {
    X<-X[idControl[cord]==1,]
  }
  #X<-rep(x,each=)

  options(na.action = currentOPTs$na.action)
  
  ### it might not be a matrix if it is just one. 
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  p<-ncol(X)
  if(!is.null(colnames(X))) namesX<-colnames(X)
  else if(p>0) namesX <- paste("var",seq(1,p),sep="")
  
  cause<-((cause==event)*1+(cause!=event & cause!=0)*2)
  
  if (Truncation){
    
    d1<-data.table::data.table(entry,exit,cause, cluster, idControl)
    data.table::setkey(d1, cluster, idControl)
    maxunexp<-d1[,max(length(idControl))-1, by=cluster][,max(V1)]
    d2<-data.table::dcast(d1, as.formula(paste("cluster","idControl", sep="~")),
                          value.var = c("entry","exit","cause"))
    if (ncol(X)!=0)  d2<-data.table::data.table(d2, X)
    data.table::setnames(d2, c("entry_1","exit_1","cause_1"), c("eentry","eexit","ecause"))
    colentry = paste("entry", (1:maxunexp)+1, sep = "_")
    colexit = paste("exit", (1:maxunexp)+1, sep = "_")
    colcause = paste("cause", (1:maxunexp)+1, sep = "_")
    if (ncol(X)==0) {
      d3<-data.table::melt(d2, id.vars=c("cluster","eentry","eexit","ecause"),
                           measure=list(colentry, colexit, colcause), value.name=c("uentry","uexit","ucause"), 
                           variable.name="idControl")
    } else d3<-data.table::melt(d2, id.vars=c("cluster","eentry","eexit","ecause", namesX),
                                measure=list(colentry, colexit, colcause), value.name=c("uentry","uexit","ucause"), 
                                variable.name="idControl")
  } else {
    
    d1<-data.table::data.table(exit,cause, cluster, idControl)
    data.table::setkey(d1, cluster, idControl)
    maxunexp<-d1[,max(length(idControl))-1, by=cluster][,max(V1)]
    d2<-data.table::dcast(d1, as.formula(paste("cluster","idControl", sep="~")),
                          value.var = c("exit","cause"))
    if (ncol(X)!=0) d2<-data.table::data.table(d2,X)
    data.table::setnames(d2, c("exit_1","cause_1"), c("eexit","ecause"))
    colexit = paste("exit", (1:maxunexp)+1, sep = "_")
    colcause = paste("cause", (1:maxunexp)+1, sep = "_")
    if (ncol(X)==0) {
      d3<-data.table::melt(d2, id.vars=c("cluster","eexit","ecause"),
                         measure=list(colexit, colcause), value.name=c("uexit","ucause"), 
                         variable.name="idControl")
    } else d3<-data.table::melt(d2, id.vars=c("cluster","eexit","ecause", namesX),
                                measure=list(colexit, colcause), value.name=c("uexit","ucause"), 
                                variable.name="idControl")
  }
  
  d3<-d3[!(is.na(uexit) & is.na(ucause))]
  d3$subj <- 1:nrow(d3)
  #strata<-d3[, get(strata), with=FALSE]

    # stacked<-prep.glm.comprisk(d, time="time", cause ="cause", times = times2, 
    #                            cens.code = 0, censmod = 0)
  #mm <- c()
  #browser()
  time.points<-matrix(time.points, ncol=1)
  #for (h in time.points) {
    if (Truncation) {
      mm<-plyr::adply(time.points,1, function(x) { 
        i2out <- prep.match.comp.risk(d3,times = x,
                                      eentrytime = "eentry" , uentrytime="uentry",
                                      eexittime = "eexit", uexittime="uexit", 
                                      ecause = "ecause", ucause="ucause",
                                      cens.code = cens.code, cens.formula = cens.formula, strata = strata)
        
        Rt <- (i2out[,eexit ] <x) * (i2out[, ecause] == 1)-(i2out[,uexit ] <x) * (i2out[, ucause] == 1)
        i2out[, ci:=ifelse(ecause!=0 & ucause!=0,Inf,ifelse(ecause==0, 
                                                            ifelse(ucause==0,min(eexit,uexit),eexit),
                                                            uexit)),by=subj]
        nocens <- (i2out[, ci] >= x) 
        return(i2out)
        },.id = NULL)
    } else {
      mm<-plyr::adply(time.points,1, function(x) { 
        #browser()
      i2out <- prep.match.comp.risk(d3,times = x,
                                    eentrytime = NULL, uentrytime=NULL,
                                    eexittime = "eexit", uexittime="uexit", 
                                    ecause = "ecause", ucause="ucause",
                                    cens.code = cens.code, cens.formula = cens.formula, strata = strata)
      Rt <- (i2out[,eexit ] <x) * (i2out[, ecause] == 1)-(i2out[,uexit ] <x) * (i2out[, ucause] == 1)
      i2out[, ci:=ifelse(ecause!=0 & ucause!=0,Inf,ifelse(ecause==0, 
                                                          ifelse(ucause==0,min(eexit,uexit),eexit),
                                                          uexit)),by=subj]
      nocens <- (i2out[, ci] >= x) 
      i2out<-cbind(i2out, Rt, h=x, nocens)
      return(i2out)
      },.id = NULL)
      
    }
    
    ### ------------- version in time ----------------------------------
    #Rt <- (i2out[,eexit ] <h) * (i2out[, ecause] == 1)-(i2out[,uexit ] <h) * (i2out[, ucause] == 1)
    
    ### ------------- keep the same value in time -----------------------
    # Rt <- (i2out[,eexit ] <h) * (i2out[, ecause] == 1)-(i2out[,uexit ] <h) * (i2out[, ucause] == 1)
    # Rtfix<- ifelse(Rt==0 &  i2out[,eexit ] <h & i2out[,uexit ] <h,
    #                ifelse(i2out[,eexit ] < i2out[,eexit ],1,-1), Rt)
    # Rt<-Rtfix
    
    
    #nocens <- ((i2out[, eexit] < h)  | (i2out[, uexit] < h))
    #mm <- rbind(mm, cbind(i2out, Rt, h, nocens))
  data.table::setDT(mm)
  mm[, clust.num:=as.numeric(as.factor(cluster))]
  mm[, weightedo:=Rt*weights]
  mm<-mm[order(clust.num,subj, h)]
  
  return(mm)
} 
###}}} compdata

###{{{ prep.match.comp.risk
prep.match.comp.risk<-function (data, times = NULL,
                                eentrytime = NULL, uentrytime=NULL,
                                eexittime = "eexit", uexittime="uexit", 
                                ecause = "ecause", ucause="ucause",
                                cname = "cweight", tname = "tweight",
                                strata = strata, 
                                nocens.data = TRUE, cens.formula = NULL, cens.code = 0, prec.factor = 100, 
                                trunc.mintau = FALSE) {
  #browser()
  #suppressMessages(require(data.table))
  if(!data.table::is.data.table(data)) data.table::setDT(data)
  out<-copy(data)
  pairexittimes<-c(eexittime, uexittime)
  #pairentrytimes<-c(eentrytime, uentrytime)
  if (is.null(times)) {
    times <- max(out[, match(pairexittimes, names(out)), with=FALSE])
  }
  if (is.null(eentrytime)) {
    if (!is.null(uentrytime)) warning("no delayed entry in exposed, assumed also for the unexposed")
    #eentrytime<-rep(0,out[,.N])
  } else {
    #warning("delayed entry in exposed, assumed also for the unexposed")
    eentrytime<-out[, pmin(get(eentrytime), get(uentrytime))]
  }
  
  mtt <- max(times)
  prec.factor <- 100
  prec <- .Machine$double.eps * prec.factor
  
  trunc.model <- cens.model <- NULL
  #browser()
  # out[, pcause:=ifelse(ecause!=0 & ucause!=0,1,0)]
  # out[, pexittime:=ifelse(pcause==0,ifelse(ecause==0,
  #                                          ifelse(ucause==0,pmin(get(eexittime), get(uexittime)),get(eexittime)),
  #                                          get(uexittime)),
  #                         pmax(get(eexittime), get(uexittime)))]
  # actual ones:
  out[, pexittime:=pmin(get(eexittime), get(uexittime))]
  out[, pcause:=ifelse(pexittime==get(eexittime),ecause,ucause)]

    #out[pexittime==0, pexittime:=pexittime+runif(1,0,0.002)] 
  if (is.null(cens.formula)) {
    if (is.null(strata)) {
      if (!is.null(eentrytime)) {
        surv.trunc <- survfit(Surv(-out[, pexittime], -eentrytime + prec, rep(1, nrow(out))) ~ 1) ##survival curve with:
        ## entry time the negative of the minimum exit time in the pair
        ## exit time the negative of the maximum entry time in the pair (+noise to avoid same values)
        ## everyone has delayed entry, so status=1 for all
        ## no covariates influences the delay entry (simple KM)
        trunc.dist <- summary(surv.trunc)
        trunc.dist$time <- rev(-trunc.dist$time) #reverse the time
        trunc.dist$surv <- c(rev(trunc.dist$surv)[-1],1) #reverse the surv, delete the first 0 and add 1 at the end
        #in this way you have like a cumulative hazard 
        if (trunc.mintau == TRUE)
          Lfit <- timereg::Cpred(cbind(trunc.dist$time, trunc.dist$surv), pmin(mtt, out[, pexittime]))
        else Lfit <- timereg::Cpred(cbind(trunc.dist$time, trunc.dist$surv),out[, pexittime]) #predict the cumulative hazard value according 
        #to the exit times of the pair: the distribution of the cumulative hazard is defined on the range of the delay entry.
        Lw <- Lfit[, 2] #these are the weights due to delay entry. 
        cformula<-"Surv(eentrytime, pexittime, pcause==0)~+1"
      }
      else {
        Lw <- 1 # truncation weights
        cformula<-"Surv(pexittime, pcause==0)~+1"
      }
      ud.cens <- survfit(as.formula(cformula), data =out)
      Gfit <- cbind(ud.cens$time, ud.cens$surv)
      Gfit <- rbind(c(0, 1), Gfit) ## Gfit is the censoring function
      Gcx <- timereg::Cpred(Gfit, pmin(mtt, out[, pexittime]), strict = TRUE)[,2] #predict the censoring function at given timepoints, get just the function value, not the time
      weights <- 1/(Lw * Gcx)
      cweights <- Gcx
      tweights <- Lw
    }
    else {
      vstrata <- as.numeric(out[, strata])
      weights <- rep(1, nrow(out))
      cweights <- rep(1, nrow(out))
      tweights <- rep(1, nrow(out))
      for (i in unique(vstrata)) {
        who <- (vstrata == i)
        if (sum(who) <= 1) 
          stop(paste("strata", i, "less than 1 observation\n"))
        outs <- subset(out, who)
        if (!is.null(eentrytime)) {
          eentrytimes <- eentrytime[who]
          surv.trunc <- survfit(Surv(-outs[, pexittime],
                                     -eentrytimes + prec, rep(1, nrow(outs))) ~ +1)
          trunc.dist <- summary(surv.trunc)
          trunc.dist$time <- rev(-trunc.dist$time)
          trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
          if (trunc.mintau == TRUE)
            Lfit <- timereg::Cpred(cbind(trunc.dist$time, trunc.dist$surv),
                          pmin(mtt, outs[, pexittime]))
          else Lfit <- timereg::Cpred(cbind(trunc.dist$time, trunc.dist$surv),
                             outs[, pexittime])
          Lw <- Lfit[, 2]
          cformula<-"Surv(eentrytimes, pexittime, pcause==0)~+1"
        }
        else {
          Lw <- 1
          cformula<-"Surv(pexittime, pcause==0)~+1"
        }
        ud.cens <- survfit(as.formula(cformula), data=outs)
        Gfit <- cbind(ud.cens$time, ud.cens$surv)
        Gfit <- rbind(c(0, 1), Gfit)
        Gcx <- timereg::Cpred(Gfit, pmin(mtt, outs[, pexittime]),strict = TRUE)[, 2]
        weights[who] <- 1/(Lw * Gcx)
        cweights[who] <- Gcx
        tweights[who] <- Lw
      }
    }
  } else {
    X <- model.matrix(cens.formula, data = out)[, -1, drop = FALSE]
    if (!is.null(eentrytime)) {
      trunc.model <- coxph(Surv(-out[, pexittime], -eentrytime + prec, rep(1, nrow(out))) ~ X)
      baseout<- basehaz(trunc.model, centered = FALSE)
      baseout <- cbind(rev(-baseout$time), rev(baseout$hazard))
      if (trunc.mintau == TRUE)
        Lfit <- timereg::Cpred(baseout, pmin(mtt, out[, pexittime]))[,-1]
      else Lfit <- timereg::Cpred(baseout, out[, pexittime])[, -1]
      RR <- exp(as.matrix(X) %*% coef(trunc.model))
      Lfit <- exp(-Lfit * RR)
      Lw <- Lfit
      cformula<-"Surv(eentrytime, pexittime, pcause==0)~+X"
    }
    else {
      Lw <- 1
      cformula<-"Surv(pexittime, pcause==0)~+X"
    }
    cens.model <- coxph(as.formula(cformula), data=out)
    baseout <- basehaz(cens.model, centered = FALSE)
    baseout <- cbind(baseout$time, baseout$hazard)
    Gfit <- timereg::Cpred(baseout, pmin(mtt, out[, pexittime]), strict = TRUE)[,2]
    RR <- exp(as.matrix(X) %*% coef(cens.model))
    Gfit <- exp(-Gfit * RR)
    weights <- 1/(Lw * Gfit)
    cweights <- Gfit
    tweights <- Lw
  }
  
  out[, (cname):=cweights]
  out[, (tname):=tweights]
  if (!is.null(eentrytime)) {
    mint <- min(tweights)
    maxt <- min(tweights)
    if (mint < 0 | mint > 1)
      warning("min(truncation weights) strange, maybe prec.factor should be different\n")
    if (maxt < 0 | maxt > 1)
      warning("max(truncation weights) strange, maybe prec.factor should be different\n")
  }
  if ("weights" %in% names(out)) {
    warning("Weights in variable 'weights_' \n")
    out[, weights_ :=weights]
  }
  else out[, weights:=weights]
  if ("cw" %in% names(out)) {
    warning("cw weights in variable 'cw_' \n")
    out[, cw_ :=1]
  }
  else out[, cw:=1]
  if (nocens.data) {
    med <- out[((pexittime > mtt & pcause== cens.code) | pcause != cens.code)]
    out<-copy(med)
  }
  attr(out, "trunc.model") <- trunc.model
  attr(out, "cens.model") <- cens.model
  return(out)
}

###}}} prep.match.comp.risk


##### to predict excess risk -------------
### risk prediction ---------------
Ft <- function(p,times,formula,newdata){
  #browser()
  lt<-length(times)
  mm<-model.matrix(as.formula(paste0("~",formula,collapse = "")), data=newdata)
  if (!all.equal(colnames(mm), names(p))) stop("check the formula, the new data and the model: you might have a different level order") 
  zp <- mm %*% p
  #if((length(times)+ncol(xpred))!=length(p)) stop("The variables in the model and the variables for prediction are different")
  #xp <- sum(xpred *  p[-(1:lt)])
  lam <- exp(zp)
  if(ncol(newdata)>1) rownames(lam)<-paste0(newdata[,1],";", newdata[,2])
  else rownames(lam)<-paste0(newdata[,1],";")
  return(lam)
}

##### new CIFpredict --- to use (new name for package)
### {{{ ecif.pred
##' Excess CIF prediction based on newdata and geepack::geese estimates.
##' @param model geese object. To be defined if coefs and vcov are null.
##' @param times vector of timepoints as the one used to estimate the GEE model
##' @param formula model formula
##' @param dataset new data
##' @param strata.levels if CIF predicted for different strata, define strata levels
##' @param coefs coefficient estimates (\code{model$beta}). To be specified if model is NULL
##' @param vcov coefficient variance and covariance matrix (\code{model$vbeta}). To be specified if model is NULL, together with coefs.
##' @import tidyr
##' @importFrom dplyr mutate
##' @examples 
##' dcif<-sim.data.MatchCR(1000,5)
##' tp<-c(0.5,1,2,5,10,15,25)
##' setdcif1<-compcomp(timereg::Event(time=FALSE,time2=time,cause=cause)~X1+X2, 
##'                    data=dcif, cluster=i, idControl=j, time.points=tp,
##'                    cens.formula=NULL, event=1)
##' exc.cif.mod1<-geepack::geese(Rt~-1+factor(h)+X1+X2,
##'                     data=setdcif1,
##'                     family="gaussian", #error distribution
##'                     mean.link = "log", #link function for Rt
##'                     corstr="independence", #correlation structure
##'                     id=clust.num, #cluster vector
##'                     weights=weights #censoring weights 
##'                     )
##' ### prediction: 
##' af<-paste0("-1+factor(h)+X1+X2") #model formula
##' newd<-data.frame(expand.grid(h=tp,X1=c(0,1),X2=c(0.8,1.5,2.5))) # newdata
##' # define the different subjects for whom the excess risk is predicted
##' strata.levels<-factor(1:6, levels=1:6,
##'                       labels =paste0(rep("X1=",6),
##'                       expand.grid(X1=c(0,1),X2=c(0.8,1.5,2.5))[,1], 
##'                       rep(", X2=",6),
##'                       expand.grid(X1=c(0,1),X2=c(0.8,1.5,2.5))[,2]))
##' pred.exc.cif<-ecif.pred(exc.cif.mod1,times = tp,dataset = newd,
##'                         formula = af, strata.levels = strata.levels)
##' head(pred.exc.cif,8)
##' @return dataset with predicted values; ready to be used with ggplot2
##' @export 
ecif.pred<-function(model=NULL,times,formula,dataset,strata.levels=NULL,
                      coefs=NULL, vcov=NULL){
  #browser()
  if( is.null(model) & is.null(coefs) & is.null(vcov)) stop("please define model or coefs+vcov")
  lt<-length(times)
  if(!is.null(model)){
    coefs<-model$beta
    vcov<-model$vbeta
  }
  
  estimate.output<-lava::estimate(model, function(p) Ft(p,times = times,
                                                        formula = formula
                                                        ,newdata = dataset), 
                                  coef=coefs, vcov = vcov
                                  #, labels = rownames(check)
  )
  
  
  if (!is.null(strata.levels)) {
    #browser()
    out.list<-estimate.output[[2]][,-5]
    strata.stop<-(1:length(strata.levels))*length(times)
    out<-do.call("rbind",
                 plyr::llply(strata.stop, function(l) {
                   #browser()
                   outl<-rbind(rep(0,4),out.list[(l-(lt-1)):l,])
                   return(outl)
                 }
                 ))
    out<-data.table(time=rep(c(0,times),length(strata.levels)),strata=rep(strata.levels, each=length(times)+1),out)
    colnames(out)[-c(1,2)]<-c("cif","se.cif","lower.ci","upper.ci")
  } else {
    out<-data.table(rbind(rep(0,5),cbind(time=times,estimate.output[[2]][,-5])))
    colnames(out)[-1]<-c("cif","se.cif","lower.ci","upper.ci")
  }
  
  return(data.table(out))
  
}

### }}} ecif.pred


### {{{ ecif.coef
##' Excess CIF prediction based on newdata and geepack::geese estimates.
##' @param model geese object. To be defined if coefs and vcov are null.
##' @param times vector of timepoints as the one used to estimate the GEE model
##' @param link link used to estimate the model. Choose between c("log","logit","identity")
##' @importFrom rlang .data
##' @examples 
##' dcif<-sim.data.MatchCR(1000,5)
##' tp<-c(0.5,1,2,5,10,15,25)
##' setdcif1<-compcomp(timereg::Event(time=FALSE,time2=time,cause=cause)~X1+X2, 
##'                    data=dcif, cluster=i, idControl=j, time.points=tp,
##'                    cens.formula=NULL, event=1)
##' exc.cif.mod1<-geepack::geese(Rt~-1+factor(h)+X1+X2,
##'                              data=setdcif1,
##'                              family="gaussian", #error distribution
##'                              mean.link = "log", #link function for Rt
##'                              corstr="independence", #correlation structure
##'                              id=clust.num, #cluster vector
##'                              weights=weights #censoring weights 
##'                              )
##' ecif.coef(exc.cif.mod1,times = tp, link = "log")
##' @return table with coefficient estimates, sandwich standard error, function of coefficient for interpretation and p-value.
##' @export
ecif.coef<-function(model, times, link="log", level=0.95){
  if (is.null(times)) stop("Vector of time points needed")
  summ<-summary(model)$mean
  p<-nrow(summ)-length(times)
  res<-summ[length(times)+(1:p), -3]
  lev<-rownames(res)
  z<-qnorm(1-(1-level)/2)
  
  f.beta<-function(beta,f.link){
    if (f.link %in% c("log","logit")) {
      f<-exp(beta)
    }
    else stop("link not supported")
    return(f)
  }

  res<-dplyr::mutate(res, 
                     lower.ci=(.data$estimate-z*.data$san.se),
                     upper.ci=(.data$estimate+z*.data$san.se))
  
  res<-dplyr::mutate(res,
                     coef=sprintf("%0.4f", .data$estimate),
                     "se(coef)"=sprintf("%0.4f", .data$san.se),
                     "f(coef)"=sprintf("%0.4f",f.beta(.data$estimate,link)),
                     lower.f.ci=sprintf("%0.4f",f.beta(.data$lower.ci,link)),
                     upper.f.ci=sprintf("%0.4f",f.beta(.data$upper.ci,link)),
                     p=sprintf("%0.4f", .data$p)) 
  rownames(res)<-lev
  res<-res[, c("coef","se(coef)","f(coef)","lower.f.ci","upper.f.ci","p")]

  return(res)
}

### }}} ecif.coef

