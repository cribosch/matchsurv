# matchsurv

## Installation 
devtools::install_github("cribosch/matchsurv")


## What can you do?
The package provides functions to estimate the cumulative excess risk for exposed individuals when matched survival data are available. 

For each exposed individual we have a defined number of unexposed individuals, matched according to some relevant factors (the number of unexposed individuals per exposed can be different). 

Here it is an example of the data
d<-data.sim(5000,5)
other options will let you chose for a competing risk setting; when nullmod=TRUE no covariates are simulated.
