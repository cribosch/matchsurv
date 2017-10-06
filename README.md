
<!-- README.md is generated from README.Rmd. Please edit that file -->
matchsurv
=========

The goal of matchsurv is to estimate the cumulative excess risk for exposed individuals when matched survival data are available.

Installation
------------

You can install matchsurv from github with:

``` r
# install.packages("devtools")
devtools::install_github("cribosch/matchsurv")
```

Example
-------

### matched survival data

``` r
library(matchsurv)
#> Loading required package: survival
#> Loading required package: timereg
#> Loading required package: mets
#> Loading required package: lava
#> lava version 1.5.1
#> mets version 1.2.2
```

For each exposed individual we have a defined number of unexposed individuals, matched according to some relevant factors (the number of unexposed individuals per exposed can be different).

Here it is an example of the data:

``` r
d<-data.sim(5000,5)
head(d,10)
#>        time status expo id j x z       cc
#> 1  63.01469      0    0  1 2 0 1 63.04353
#> 2  64.42602      0    0  2 2 0 1 65.81813
#> 3  65.95442      0    0  3 2 0 1 58.05087
#> 4  32.38016      1    0  4 2 0 0 65.50442
#> 5  79.61409      0    0  5 2 1 1 63.75141
#> 6  73.88104      0    0  6 2 1 1 70.98413
#> 7  75.11966      0    0  7 2 1 0 67.36123
#> 8  62.00047      0    0  8 2 0 1 66.53394
#> 9  16.02534      1    0  9 2 1 1 71.87354
#> 10 73.51231      0    0 10 2 1 0 61.99114
```

`matchsurv::data.sim` let you simulate some matched survival data; `competing=TRUE` will let you chose for a competing risk setting; when `nullmod=TRUE` no covariates are simulated.

This is a basic example which shows you how to:

1.  estimate the proportional excess model
2.  visualize your results - coefficient estimates and cumulative excess hazard plot

### model estimate

First you need to set up your data in order to estimate the model

``` r
example("compdata")
```

`strata()` are not needed in here; you'll choose after how to specify your model. It is a good idea to use all the possible covariates at this step. New variables will be created, you'll need them in the following step.

Then you can estimate your model:

``` r
example("matchpropexc")
```

model results
-------------

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
