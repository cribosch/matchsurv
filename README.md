
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
#> 1  65.05970      0    0  1 2 0 1 65.90091
#> 2  79.55744      0    0  2 2 1 0 65.80985
#> 3  75.29401      0    0  3 2 1 0 70.66296
#> 4  77.14508      0    0  4 2 1 1 58.69530
#> 5  75.18099      0    0  5 2 1 1 65.31134
#> 6  73.23426      0    0  6 2 1 1 65.37848
#> 7  67.10663      0    0  7 2 0 1 67.99293
#> 8  79.14640      0    0  8 2 1 0 69.37830
#> 9  69.63919      0    0  9 2 0 0 67.98359
#> 10 74.30575      0    0 10 2 1 0 67.78729
```

Or when working with cumulative incidence functions:

``` r
dcif<-sim.data.MatchCR(1000,5)
head(dcif, 10)
#>      i j expo X1        X2      agee     entry      exit      time cause
#>  1:  1 1    1  0 1.3645013 14.568374 14.568374 29.034426 14.466052     1
#>  2:  2 1    1  1 1.5996915 18.028515 18.028515 24.116528  6.088014     1
#>  3:  3 1    1  1 1.3139423 17.485810 17.485810 18.638367  1.152557     2
#>  4:  4 1    1  1 1.7093352  7.966872  7.966872  9.800529  1.833657     1
#>  5:  5 1    1  1 1.4271260 24.179491 24.179491 50.000000 25.820509     0
#>  6:  6 1    1  0 1.3220371 18.428603 18.428603 39.741398 21.312795     1
#>  7:  7 1    1  0 1.4367239  5.906687  5.906687  6.382655  0.475968     1
#>  8:  8 1    1  0 1.1640960 16.500282 16.500282 28.737923 12.237641     1
#>  9:  9 1    1  1 1.8546899 14.050822 14.050822 21.504109  7.453287     1
#> 10: 10 1    1  1 0.9732946 10.207745 10.207745 21.072370 10.864625     1
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

### Excess cumulative incidence

To visualize the estimated effects:

``` r
example("ecif.coef")
```

To predict the excess cumulative incidence for different covariate values:

``` r
example("ecif.pred")
```

Further info on plotting the results available in `vignette("matchsurv")`
