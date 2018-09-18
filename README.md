
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
#> Loading required package: lava
#> lava version 1.6.3
#> Loading required package: timereg
#> Loading required package: survival
#> Loading required package: mets
#> mets version 1.2.4
#> Loading required package: data.table
#> Loading required package: plyr
```

For each exposed individual we have a defined number of unexposed individuals, matched according to some relevant factors (the number of unexposed individuals per exposed can be different).

Here it is an example of the data:

``` r
d<-data.sim(5000,5)
head(d,10)
#>        time status expo id j x z       cc
#> 1  77.95430      0    0  1 2 1 0 74.04276
#> 2  66.99011      0    0  2 2 0 1 64.23351
#> 3  61.16902      0    0  3 2 0 1 67.23846
#> 4  37.74204      1    0  4 2 1 1 69.15333
#> 5  66.56251      0    0  5 2 0 1 72.16574
#> 6  74.32605      0    0  6 2 1 1 67.13336
#> 7  66.52239      0    0  7 2 0 1 71.54178
#> 8  65.62297      0    0  8 2 0 1 58.51820
#> 9  74.34668      0    0  9 2 1 0 62.66389
#> 10 29.54146      1    0 10 2 1 1 72.88285
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
