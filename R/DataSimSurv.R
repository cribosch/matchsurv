## This file contains the function to simulate the data.
## if option Trunc you can decide if to simulate left-truncated data or not.

##maybe you can add it later to the package. Let it here for now and recall it into your simulation file. 

### {{{ dataset.sim

dataset.sim<-function(alpha0,lambda0,n,k,gamma,betacoef,
                      cens=FALSE, tcens, Trunc=FALSE){
    Z<-rbinom(n,1,0.5)
    X1<-rbinom(n,1,0.5)
    X2 <- rbinom(n,1,0.5)
    
    Zk<-rep(Z, each=k)
    X1k<-rep(X1, each=k)
    X2k<-rep(X2, each=k)
    
    #controls
    lamk <-alpha0*exp(gamma*Zk) 
    t_star_k <- rexp(n*k)/lamk
    if (Trunc) t_start_k <- runif(n*k,0,3)
    ck <- runif(n*k,0,12)
    
    #cases 
    lamc <- alpha0*exp(gamma*Z)+lambda0*exp(betacoef[1]*X1+betacoef[2]*X2)


    t_star_c <- rexp(n)/lamc
    if (Trunc) t_start_c <- runif(n,0,3)
    cc <- runif(n,0,12)

    if (Trunc)t_start <- c(t_start_k,t_start_c)
    t_star <-c(t_star_k,t_star_c)
    t_cens <-c(ck,cc)
    i <- c(rep(seq(1:n),each=k),1:n)
    j <- c(rep(2:(k+1),n), rep(1,n))
    Z <- c(Zk, Z)
    X1<- c(X1k, X1)
    X2<- c(X2k, X2)

    time <- pmin(t_star,t_cens)
    status<- (time==t_star)*1

    if (Trunc) {
        t_end <- t_start+time
    }
    else t_end <- time

    if (cens==TRUE){
        status <- ifelse(t_end<=tcens,status,0)
        t_end<- ifelse(t_end<=tcens, t_end, tcens)
    }

    if (!Trunc) t_start=rep(0,length(time))

    dataframe <- data.frame(i=i,
                            j=j,
                            X1=X1,
                            X2=X2,
                            Z=Z,
                            start=t_start,
                            time=time,
                            end=t_end,
                            status=status)
    return(dataframe)
                          
}


### This is another function to simulate data based on pc.hazard function of timereg package. It is possible to simulate
### also cometing risk data. It is not possible to define the hazard. 

####{{{ sim.data.MatchH
##' Function to simulate macthed survival data
##' @param nca number of exposed individuals
##' @param ncont number of unexposed individuals for eac exposed (fixed number)
##' @param competing if TRUE a competing cause is considered
##' @param nullmod if TRUE no covariates are simulated. By default 3 variables are simulated: two binomals and one continuous. 
##' @author Cristina Boschini
##' @export
 sim.data.MatchH<- function(nca, #number of cases
                       ncont, #number of controls
                       competing=FALSE, #with a competing event
                       nullmod=FALSE){
  #browser()

  n <- nca
  
  ### covariates
  #x <- rbinom(n,1,0.6)
  #sa <- 7
  entry <- alder <- runif(n)*20+5
  #other covariates
  if(!nullmod) {
    z <- rbinom(n,1,0.8)
    cc <- rnorm(n,65,4)
  }
  
  ## background hazards
  haz1 <- 0.07/40  #first 50 years: 0.023, other 4 years is 0.154
  haz2 <- 0.3/80
  ## back 1
  back1 <- rbind(c(0,0),c(40,40*haz1), c(80,haz1*40+40*haz2))
  ## back 2
  if (competing)  back2 <- rbind(c(0,0),c(40,0.25), c(80,0.35)) #hazard of the other cause
  
  # time/status for cause 1
  cont1 <- lapply(1:ncont, function(y) {
    data.frame(timereg::pc.hazard(back1,n=n,entry=entry),  #time for unexposed with background=1
               j=y+1,
               id=1:n,
               expo=0)
  }
  )
  
  cont1 <- do.call("rbind", cont1)
  cont1$codur= cont1$time-cont1$entry
  cont1$costatus <- cont1$status
  
  
  if (competing) {
    
    contcomp <- lapply(1:ncont, function(y) {
      data.frame(timereg::pc.hazard(back2, n=n, entry=entry))
    }
    )
    
    contcomp <- do.call("rbind", contcomp)
    
    agee <- pmin(cont1$time, contcomp$time) 
    cont1$costatus <- ifelse(cont1$time<contcomp$time, cont1$status, 2*contcomp$status)
    cont1$codur <- agee-entry
    cont1$agestop <- agee
  }
  
  ## exposed
  excess1 <- rbind(c(0,0),c(5,0.20), c(40,0.30))
  case11 <- timereg::pc.hazard(back1,n=n,entry=entry)
  if(!nullmod){
    rr <- exp(z*0.5)
    case12 <- timereg::pc.hazard(excess1,rr=rr)
  } else case12 <- timereg::pc.hazard(excess1,n=n)
  
  case11$dur <- case11$time-case11$entry
  case11$dur  <- pmin(case11$dur,case12$time) 
  case11$status <- ifelse(case11$dur<case12$time,case11$status,case12$status)
  case11$id <- 1:n
  case11$j <- 1
  case11$expo <- 1
  
  if (competing) {
    excess2 <-t(c(1,0.8)*t(excess1))
    if(!nullmod) {
      rr <- exp(z*0.25)
      case22 <- timereg::pc.hazard(excess2,rr=rr)
    } else  case22 <- timereg::pc.hazard(excess2,n=n)
    case11$status <- ifelse(case11$dur<case22$time, case11$status, 2*case22$status)
    case11$dur <- pmin(case11$dur, case22$time)
  }
  
  ## complete dataset
  dd <- data.frame(time=c(cont1$codur,case11$dur),
                   status=c(cont1$costatus, case11$status),
                   expo=c(cont1$expo, case11$expo),
                   id=c(cont1$id, case11$id),
                   j=c(cont1$j, case11$j)
  )
  if (!nullmod) {
    dd <- data.frame(dd,
                     #x=c(rep(x,ncont+1)),
                     z=c(rep(z,ncont+1)),
                     cc=c(rep(cc,ncont+1))
    )
  }
  return(dd)
 } 
 ###} sim.data.MatchH
 ###} dataset.sim

