## This file contains the function to simulate 
## matched data in a competing risk setting, based
## on the cumulative incidence scale

### global variables
if(getRversion() >= "2.15.1")  utils::globalVariables(c("exit","i","j","agee"))  

####{{{ simMatchCR
##' Function to simulate matchced survival data in competing risk setting based on cumulative incidence functions
##' @param nca number of exposed individuals
##' @param ncont number of unexposed individuals for each exposed (fixed number; by default ncont=5)
##' @param cifs list argument. Starting CIFs for: 1) 1st cause background risk, 2) 2nd cause background risk, 3) 1st cause excess risk, 4) 2nd cause excess risk. If cifs=NULL default quantities are used. 
##' @param gammax coefficients for the excess risk model. Max 2 values to be specified (the first covariate is by default binomial, the second is log normal), Default values: \eqn{\gamma_1=0.1, \; \gamma_2=-0.2}. Note: some value might not work.
##' @param mean.link link function for the gee model. The data are simulating accordingly to the link that will be used to estimate the excess cif model. Two values supported "id" or "log".
##' @param bias to simulate data with excess risk factor correlated with age
##' @param cens TRUE when simulating data with uniform censoring time (it's the same for exposed and unexposed in the same cluster to mimic the registry structure) 
##' @param cens.dist TRUE when censoring time depends on X1. FALSE if no censoring or independent censoring times
##' @param age.expo NULL if age at entry randomly generated from uniform(0,20). If fix age, specify single number (it can be 0, no delay entry)
##' @param print.cifs if information about the given cifs is needed
##' @importFrom utils tail
##' @importFrom stats model.matrix as.formula
##' @author Cristina Boschini
##' @return The function returns by default a dataset. if print.cifs=TRUE, the function returns a list withe the simulated dataset and a list of the startingc cumulative incidence functions.
##' @export

sim.data.MatchCR<-function(nca,
                           ncont=5,
                           cifs=NULL,
                           gammax=c(0.1,-0.2), #by default simulate two covariates, z can be 0,1,2. 
                           mean.link="log", 
                           bias=FALSE,
                           cens=FALSE,
                           cens.dist=FALSE,
                           age.expo=NULL, #age at entry for exposed if null, uniform (0,20). specify value otherwise. it can be 0
                           print.cifs=FALSE
){
  #browser()
  if (is.null(cifs)) {
    ecif1<-cbind(c(0,0.2,0.5,1:30),
                 c(0,0.074,0.098,0.117,0.117,
                   0.124,0.136,0.138,0.135,0.144,0.148,0.15,0.156,
                   0.171,0.173,0.161,0.161,0.157,0.157,0.157,0.157,
                   0.157,0.15,0.15,0.15,0.15,0.15,0.148,0.148,0.148,0.148,0.148,0.148))    
    ecif2<-cbind(c(0,0.2,0.5,1:30), c(rep(0,2),rep(0.01,2), 0.02, rep(0.03,11),rep(0.04,17)))
    ucif1<-cbind(0:50,c(0,0.055,0.095,0.127,0.166,0.196,
                        0.23,0.259,0.279,0.309,0.333,0.36,
                        0.387,0.417,0.444,0.479,0.507,0.529,
                        0.545,0.565,0.585,0.604,0.621,0.637,
                        0.652,0.664,0.677,0.69,0.704,0.714,
                        0.724,0.737,0.746,0.757,0.765,0.774,
                        0.784,0.793,0.805,0.813,0.824,0.832,
                        0.843,0.851,0.858,0.861,0.866,0.871,
                        0.881,0.889,0.897))
    ucif2<-cbind(0:50,c(rep(0,8),rep(0.001,2), rep(0.002, 3),
                        rep(0.003,4), rep(0.004,8), rep(0.005,4), rep(0.006,9),
                        rep(0.007,2), rep(0.008,2), rep(0.009,4), rep(0.010,5)))
    uF<-list(uF1=ucif1,uF2=ucif2)
    excess<-list(excess1=ecif1,excess2=ecif2)
    cifs<-list(uF[[1]], uF[[2]], excess[[1]], excess[[2]])
  } else {
    if (!is.list(cifs)) warning("cifs is a list with 4 arguments in this order: background risk 1st cause,
                                background risk 2nd cause, excess risk 1st cause, excess risk 2nd cause")
    uF<-list(uF1=cifs[[1]], uF2=cifs[[2]])
    excess<-list(excess1=cifs[[3]], excess2=cifs[[4]])
  } 
  
  if (bias) {
    uF1.f<-function(t){
      uF1<-(t<=10)*(t*0.01)+(t>10 & t<=25)*(0.1+(t-10)*0.03)+(t>25)*(0.55+(t-25)*0.01)
      return(cbind(time=t,cuminc=uF1))
    }
    uF1<-uF1.f(uF[[1]][,1])
    uF[[1]]<-uF1
  }
  
  
  if (!is.null(gammax)){
    if (length(gammax)>1) namesX<-c("X1","X2")
    else namesX<-"X1"
  } else if (bias) namesX<-"X1"
  
  #### exposed ------------
  expo<-plyr::rdply(nca, {
    F1<-simexpo(F0=uF,
                Fe=excess,
                model=mean.link,
                gamma = gammax,
                biass=bias,
                age.entry = age.expo)
    #browser()
    entry<-F1$entry
    if(!is.null(gammax)) {
      X<-matrix(F1$X, ncol=length(gammax))
      colnames(X)<-namesX
    }
    if(bias) {
      X<-matrix(F1$X, ncol=1)
      colnames(X)<-namesX
    }
    F1.cuminc<-F1[[1]]
    ptot<-sum(plyr::ldply(F1.cuminc, function(x) {
      tail(x[,2],1)
    })[,2])
    if (ptot<=1) {
      expotc<-data.frame(timecause(F1.cuminc[[1]],F1.cuminc[[2]],n=1,entry=NULL,ptot=ptot))
      if (cens){
        if (!cens.dist){
        censd<-rexp(1,0.02)
        } else {
          if (!("X1" %in% namesX)) stop("X1 is needed")
          censd<-rexp(1,0.01)*(X[,1]==1)+rexp(1,0.1)*(X[,1]!=1)
        }
        time<-pmin(expotc$time,censd)
        expotc$cause <-ifelse(expotc$time==time,expotc$cause,0)
        expotc$time<-time
        expotc$censd<-censd
      }
      
      e<-data.frame(expotc, j=1,expo=1, entry=entry)
      if(!is.null(gammax)) e<-cbind(e,X)
      else if (bias)  e<-cbind(e,X)
      return(e)
    }
  })
  
  
  #browser()
  expo$exit<-expo$time+expo$entry
  data.table::setDT(expo)
  if (expo[,.N]<nca) stop("ptot >1") #check2
  names(expo)[1]<-"i"
  
  #### unexposed -----------------
  unexpo <- plyr::ldply(1:ncont, function(y){
    #browser()
    entry<-expo$entry
    if(!is.null(gammax) | bias) X<-expo[, namesX, with=FALSE]
    unexpotc<-data.frame(timecause(uF[[1]],uF[[2]],n=nca,entry=entry))
    if(cens){
      censa<-expo$censd+entry
      # unexpotc$cause <-ifelse(unexpotc$time<=censa,unexpotc$cause,0)
      # unexpotc$time[unexpotc$cause==0]<-censa[unexpotc$cause==0]
      # unexpotc$censd<-expo$censd
      time<-pmin(unexpotc$time,censa)
      unexpotc$cause <-ifelse(unexpotc$time==time,unexpotc$cause,0)
      unexpotc$time<-time
      unexpotc$censd<-expo$censd
      
    }
    u<-data.frame(unexpotc,
                  entry=entry,
                  expo=0,
                  i=1:nca,
                  j=y+1)
    if(!is.null(gammax) | bias) u<-cbind(u,X)
    return(u)
  })
  
  unexpo$exit<-unexpo$time
  unexpo$time<-unexpo$exit-unexpo$entry
  data.table::setDT(unexpo)
  
  #### together ----------
  dd <- rbind(expo, unexpo)
  
  ### check list ---------
  if (dd[,sum(time<=0)]) message("Negative times simulated")
  if (dd[,sum(time<=0)]) print(dd[time<=0])
  if (dd[,sum(duplicated(cbind(i,j)))>0])  message("Not unique id")
  
  #browser()
  ### covariates ------
  dd[, agee:=entry[expo==1],by=i]
  if (!is.null(gammax)) {
    if (length(gammax)>1) data.table::setcolorder(dd, c("i", "j", "expo","X1","X2","agee","entry","exit","time","cause"))
    else data.table::setcolorder(dd, c("i", "j", "expo","X1","agee","entry","exit","time","cause"))
  } else data.table::setcolorder(dd, c("i", "j", "expo","agee","entry","exit","time","cause"))
  
  if (!print.cifs) return(dd)
  else return(list(data=dd, cifs=cifs))
}
####}}} simMatchCR

Fexpo<-function(Funexpo,
                Fexcess,
                coefs,
                entryt,
                truncprob,
                cov,
                mean.link) {
  #browser()
  if (is.list(Funexpo) & length(Funexpo)==1) {
    Funexpo<-Funexpo[[1]]
    Fexcess<-Fexcess[[1]]
  }
  ucifentry<-timereg::subdist(Funexpo, entryt)[,2]
  F1entrya<-Funexpo
  F1entrya[,1]<-F1entrya[,1]-entryt
  F1entrya[,2]<-(F1entrya[,2]-ucifentry)/truncprob
  F1entryaL<-rbind(c(0,0), F1entrya[F1entrya[,1]>0,])
  F1times<-F1entryaL[,1,drop=FALSE]
  exctimes<-Fexcess[,1, drop=FALSE]
  maxd<-tail(exctimes,1)[1]
  FeZ<-cbind(exctimes,expand.grid(Fexcess[,2],sum(coefs*cov)))
  if (mean.link=="id") FeZ<-cbind(FeZ[,1],FeZ[,2]+FeZ[,3])
  else FeZ<-cbind(FeZ[,1],FeZ[,2]*exp(FeZ[,3]))
  if( FeZ[1,2]!=0) FeZ[1,2]<-0
  #F1t<-cbind(time=F1times,Fe=F1entryaL[,2]+timereg::Cpred(FeZ,F1times,strict=FALSE)[,2])
  F1t<-cbind(time=F1times,Fe=F1entryaL[,2]+Hmisc::approxExtrap(FeZ[,1],FeZ[,2],xout = F1times)$y)
  F1t<-F1t[F1t[,1]<=maxd,]
  if (any(F1t[,2]>1)) stop("gamma too big") #check3
  #browser()
  if (any(diff(F1t[,2])<0)) stop (paste0("CIF in exposed is decreasing;entryt:",entryt)) #if any TRUE stop!!!!!! #check4
  return(F1t)
}

##### simexpo ----------------------------------------
### give F1 on the age time scale
### give Fe on the duration time scale
### Z can be continuous, default is null
### !is.null(Z) -> remember the coefficient, gamma
### link to specify in order to simulate data correctly with the model we use to estimate later 
## there is a file called simexpo that contains the same function (version 22/03)+comments with the previous versions in old folder
## they will be one binomial and the other log-normal 

simexpo<-function(F0,
                  Fe,
                  model="log",
                  gamma=NULL, 
                  #beta=FALSE, ### matching factor effect assume none
                  biass=FALSE,
                  age.entry=NULL ## when age at entry is uniform (0,1). specify value if fix age.
) {
  #browser()
  
  if (!is.list(F0)) {
    F0<-list(F0=F0)
    Fe<-list(Fe=Fe)
  }
  
  twocauses<-TRUE
  if (length(F0)==1) twocauses<-FALSE
  
  twocovs<-FALSE
  if(!is.null(gamma)) {
    X1<-rbinom(1,1,0.5)  # Z<-1 #
    if (length(gamma)>1) {
      X2<-rlnorm(1,mean=0.3,sd=0.2)
      twocovs<-TRUE
    }
  } else if (biass) X1<-rbinom(1,1,0.5)
  
  # matching factor
  # if(!is.null(beta)){
  #   Z<-rbinom(1,1,0.35)
  # }
  
  if (biass){ 
    age<-(runif(1)*2)*(X1==1)+(X1==0)*((15+runif(1)*3))
  } else if (is.null(age.entry)) {
    age<-runif(1,0,20)
  } else age<-age.entry
  
  
  entry<-age
  
  #browser()
  
  # no matchinf factor effect assumed
  # if(!is.null(beta)){
  #   ti<-F0[[1]][,1]
  #   F0X1<-cbind(ti,expand.grid(F0[[1]][,2],sum(beta*X))) ## matching factor affects only first xcause
  #   F0X1<-cbind(F0X1[,1],F0X1[,2]*exp(F0X1[,3]))
  #   if( F0X1[1,2]!=0) F0X1[1,2]<-0
  #   F0[[1]]<-F0X1
  # }
  
  ptruncc<-plyr::ldply(F0, function(x) {
    ucifentry<-timereg::subdist(x, entry)[,2]
  }
  )
  ptrunc<-1-sum(ptruncc[,2])
  if (any(ptrunc>1)) stop("ptrunc >1") #check1
  
  if (!is.null(gamma)) Xv<-X1
  if (twocovs) Xv<-c(X1,X2) 
  
  if(is.null(gamma)) {
    gamma<-0
    Xv<-NULL
    if(biass) Xv<-X1
  }
  
  if (twocauses) {
    F1<-plyr::mlply(cbind(x=F0, y=Fe, g=list(gamma,0)),
                    function(x,y,g) Fexpo(x,y,g, entryt = entry, 
                                          truncprob = ptrunc, 
                                          cov=Xv, mean.link = model))
  } else F1<-Fexpo(F0,Fe,gamma, entryt = entry, 
                   truncprob = ptrunc, cov=Xv, 
                   mean.link = model)
  #browser()
  return(list(F1=F1, entry=entry, X=Xv))
}

#### timecause --------------
## this functions simulates times and causes for the two causes at the same time 

timecause<-function(Fc1,Fc2,n,entry=NULL,ptot=NULL){
  #browser()
  u<-runif(n)
  if (!is.null(entry)){
    F1entry <- timereg::subdist(Fc1,entry)[,2]
    F2entry <- timereg::subdist(Fc2,entry)[,2]
    ptrunc <- 1-F1entry-F2entry
    Fi1 <- timereg::invsubdist(Fc1,u,cond=1,entry=entry,ptrunc=ptrunc)
    Fi2 <- timereg::invsubdist(Fc2,u,cond=1,entry=entry,ptrunc=ptrunc)
  } else {
    ptrunc<-rep(1,n)
    F1entry<-rep(0,n)
    F2entry<-rep(0,n)
    Fi1<-timereg::invsubdist(Fc1,u,cond=1)
    Fi2<-timereg::invsubdist(Fc2,u,cond=1)
  }
  
  if(is.null(ptot)) ptot<-(tail(Fc1[,2],1)+tail(Fc2[,2],1)-F1entry-F2entry)/(ptrunc)
  rt <- rbinom(n,1,ptot)
  p1 <- ((tail(Fc1[,2],1)-F1entry)/ptrunc)
  rb <- rbinom(n,1,p1/ptot)
  cause <- ifelse(rb==1,1,2)
  time <- ifelse(cause==1,Fi1$time,Fi2$time)
  if (any(time>50)) print(cbind(u,cause,time)[time>50])
  cause <- rt*cause
  time[cause==0] <- min(tail(Fc1[,1],1), tail(Fc2[,1],1))
  return(cbind(time=time,cause=cause))
}


