## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

requireNamespace("devtools")

## ---- echo=TRUE, results="hide",eval=FALSE-------------------------------
#  # install.packages("devtools")
#  devtools::install_github("cribosch/matchsurv",build_vignettes = TRUE )

## ------------------------------------------------------------------------
library(matchsurv)

#other useful packages 
library(timereg)
library(geepack)
library(ggplot2) # for the plot

## ------------------------------------------------------------------------
haz.data<-sim.data.MatchH(5000,5)
head(haz.data)

## ------------------------------------------------------------------------
set.hazd<-compdata(Surv(time, status)~z+cc, clust=id, idControl=j, data=haz.data)
head(set.hazd)

## ------------------------------------------------------------------------

exc.haz.mod<-matchpropexc(Surv(entry,exit,status)~strata(z)+cc, data=set.hazd)


## ------------------------------------------------------------------------

summary(exc.haz.mod)


## ------------------------------------------------------------------------
exccumhaz(exc.haz.mod, time = seq(0,25,5))

## ------------------------------------------------------------------------
excplot(exc.haz.mod)

## ------------------------------------------------------------------------

cif.data<-sim.data.MatchCR(1000,5)
head(cif.data)

## ------------------------------------------------------------------------
tp<-c(0.5,1,2,5,10,15,25) # define vector of time points
set.cifd<-compcomp(Event(time=FALSE,time2=time,cause=cause)~X1+X2,
                   data=cif.data,
                   cluster=i, #cluster ID
                   idControl=j, #subject ID
                   time.points = tp, #time points
                   cens.formula = NULL, #censoring weight model
                   event=1 #event of interest
)
head(set.cifd)

## ------------------------------------------------------------------------
exc.cif.mod<-geese(Rt~-1+factor(h)+X1+X2,
                data=set.cifd,
                family="gaussian", #error distribution
                mean.link = "log", #link function for Rt
                corstr="independence", #correlation structure
                id=clust.num, #cluster vector
                weights=weights #censoring weights
                )


## ------------------------------------------------------------------------

af<-paste0("-1+factor(h)+X1+X2") #model formula
newd<-data.frame(expand.grid(h=tp,X1=c(0,1),X2=c(0.8,1.5,2.5))) # newdata
# define the different subjects for whom the excess risk is predicted
strata.levels<-factor(1:6, levels=1:6, 
                      labels =paste0(rep("X1=",6),
                                     expand.grid(X1=c(0,1),X2=c(0.8,1.5,2.5))[,1],
                                     rep(", X2=",6),
                                     expand.grid(X1=c(0,1),X2=c(0.8,1.5,2.5))[,2]))

pred.exc.cif<-ecif.pred(exc.cif.mod,times = tp,dataset = newd, 
                     formula = af, strata.levels = strata.levels)
head(pred.exc.cif,8)


## ------------------------------------------------------------------------
p<-ggplot()
plot.exc.cif<-p+
  geom_step(aes(x=time, y=lower.ci, alpha=0.1), data=pred.exc.cif,lty="dotted", size=0.7)+
  geom_step(aes(x=time, y=upper.ci, alpha=0.1), data=pred.exc.cif,lty="dotted", size=0.7)+
  geom_step(aes(x=time, y=cif),col="black", data=pred.exc.cif, size=0.8)+
  scale_x_continuous(name="Time since entry",limits = c(0,30),breaks = seq(0,30,5), labels=seq(5,35,5))+
  scale_y_continuous(name="Excess risk")+ 
  facet_wrap(~strata,ncol = 2)+
  geom_abline(slope = 0,intercept = 0, size=0.2)+
  theme(legend.position = 'none');plot.exc.cif

## ------------------------------------------------------------------------

ecif.coef(exc.cif.mod,times = tp, link = "log")


