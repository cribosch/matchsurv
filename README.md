
<!-- README.md is generated from README.Rmd. Please edit that file -->
matchsurv
=========

The goal of matchsurv is to estimate the cumulative excess risk for exposed individuals when matched survival data are available. Excess risk regression models are available on both hazard and cumulative incidence scale. More info in `vignette("matchsurv")`.

Installation
------------

You can install matchsurv from github with:

``` r
# install.packages("devtools")
devtools::install_github("cribosch/matchsurv")
```

and call the package with:

``` r
library(matchsurv)
#> Loading required package: survival
#> Loading required package: data.table

## other useful packages
library(timereg)
library(geepack)
library(ggplot2)
```

The package works with matched cohort data where the outcome on study is time-to-event. For each exposed individual we have a defined number of unexposed individuals, matched according to some relevant factors (the number of unexposed individuals per exposed can be different). To simulate data:

``` r
dhaz<-sim.data.MatchH(5000,5)
head(dhaz,10)
#>        time status expo id j x z       cc
#> 1  64.81160      0    0  1 2 0 1 65.23130
#> 2  74.55864      0    0  2 2 1 1 61.74196
#> 3  74.64357      0    0  3 2 1 1 67.15914
#> 4  76.72766      0    0  4 2 1 1 71.69239
#> 5  75.37207      0    0  5 2 1 1 64.30109
#> 6  74.28333      0    0  6 2 1 1 63.65761
#> 7  68.01524      0    0  7 2 0 1 65.64770
#> 8  73.63796      0    0  8 2 1 1 62.48640
#> 9  76.34294      0    0  9 2 1 1 71.05816
#> 10 76.53048      0    0 10 2 1 1 70.30413
```

Or when working with cumulative incidence functions:

``` r
dcif<-sim.data.MatchCR(1000,5)
head(dcif, 10)
#>      i j expo X1        X2      agee     entry      exit       time cause
#>  1:  1 1    1  1 1.0629432 14.415520 14.415520 24.117773  9.7022532     1
#>  2:  2 1    1  1 1.2700912 21.713826 21.713826 26.328167  4.6143416     1
#>  3:  3 1    1  0 0.9100611  9.635260  9.635260 18.244576  8.6093155     1
#>  4:  4 1    1  0 1.0178028 18.391567 18.391567 21.780796  3.3892287     2
#>  5:  5 1    1  0 1.2989273 24.295612 24.295612 33.518705  9.2230935     1
#>  6:  6 1    1  0 0.9367467 16.151441 16.151441 16.233732  0.0822906     1
#>  7:  7 1    1  0 1.4190000  7.315778  7.315778  9.150042  1.8342635     1
#>  8:  8 1    1  0 1.5175804  5.558891  5.558891 14.776881  9.2179906     1
#>  9:  9 1    1  0 1.1422343 15.021235 15.021235 15.133945  0.1127099     1
#> 10: 10 1    1  1 0.9868964  8.252132  8.252132 23.590356 15.3382242     1
```

Further info about the function options are available through `vignette()` or `example()`.

This is a basic example which shows you how to:

1.  estimate the excess risk
2.  visualize your results - coefficient estimates and cumulative excess

in both the two settings.

Model estimate
--------------

### Excess hazard

Data set up

``` r
example("compdata")
```

Model estimate

``` r
example("matchpropexc")
```

### Excess cumulative incidence

Data set up

``` r
example("compcomp")
#> Warning in example("compcomp"): 'compcomp' has a help file but no examples
```

Model estimate based on GEE function `geepack::geese()`

``` r
tp<-c(0.5,1,2,5,10,15,25)
setdcif1<-compcomp(Event(time=FALSE,time2=time,cause=cause)~X1+X2,
data=dcif, cluster=i, idControl=j, time.points=tp, cens.formula=NULL, event=1)

exc.cif.mod1<-geese(Rt~-1+factor(h)+X1+X2,
                    data=setdcif1,
                    family="gaussian", #error distribution
                    mean.link = "log", #link function for Rt
                    corstr="independence", #correlation structure
                    id=clust.num, #cluster vector
                    weights=weights #censoring weights
)
```

Results
-------

### Excess hazard

To visualize the coefficient estimates: `summary(model)`

To estimate the cumulative baseline excess hazard:

``` r
example("exccumhaz")
```

To plot the cumulative baseline excess hazard:

*Note: if your model has strata, you can chose which strata to plot (option: `stratas=`, followed by the number of the strata, the first one is number 0). You can also decide to show the relative survival (option: `relsurv=TRUE`).*

``` r
example("excplot")
```
