% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchcif.R
\name{compcomp}
\alias{compcomp}
\title{Data structured for glm approach in competing risk setting}
\usage{
compcomp(formula, data, cluster, idControl, strata = NULL, time.points,
  cens.formula = NULL, cens.code = 0)
}
\arguments{
\item{formula}{formula with 'Event' outcome (see \code{timereg} package); time stands for the start time, while time2 stands for the stop time. cause=1 will be considered as the event of interest}

\item{data}{data frame}

\item{cluster}{vector cluster indicator (one cluster for each exposed individual)}

\item{idControl}{vector control indicator (idControl==1 indicates exposed individual in cluster i)}

\item{strata}{}

\item{time.points}{vector of time points where the glm will be estimated (10 usually is a sufficient number; the more time points, the slower the glm function)}

\item{cens.formula}{useful to estimate the weights when censoring is present. no quotes, add something like ~age+year}
}
\value{
A setup dataset, ready for \code{geese}
}
\description{
Data structured for glm approach in competing risk setting
}
\author{
Cristina Boschini
}