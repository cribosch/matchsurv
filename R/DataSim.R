## This file contains the function to simulate the data.
## if option Trunc you can decide if to simulate left-truncated data or not.

##maybe you can add it later to the package. Let it here for now and recall it into your simulation file. 

##Remember: probands=cases(1); controls goes from 2 to k+1. 

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

###} dataset.sim

