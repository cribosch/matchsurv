% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchf.R
\name{predict.matchpropexc}
\alias{predict.matchpropexc}
\title{Predict cumulative excess hazard or relative survival}
\usage{
\method{predict}{matchpropexc}(object, relsurv = FALSE,
  time = object$exit, X = object$X, strata = object$strata, ...)
}
\arguments{
\item{object}{model estimated with matchpropexc}

\item{relsurv}{set to FALSE. If TRUE relative survival is returned instead of the cumulative baseline hazard}

\item{time}{specify at which time to compute the estimates}

\item{X}{define specific values for the excess covariates. By default predictworks on the original data.}

\item{strata}{default strata if the model has them.}

\item{...}{Additional arguments to lower level funtions}
}
\description{
Predict cumulative excess hazard or relative survival
}
\author{
Cristina Boschini
}
