% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.spatial.R
\name{emp.cov}
\alias{emp.cov}
\title{Collects emperical covariances, and cross covarianes}
\usage{
emp.cov(loc, Y, Bins = NULL, nBins = 10)
}
\arguments{
\item{loc}{A list of locations      (d x 2)}

\item{Y}{A list of observations   (m x 1-2)}

\item{Bins}{right locations of the bins {0},(0,Bins[1]), ... (Bins[n-1], Bins[n]), (Bins[n], \eqn{\inf})}

\item{nBins}{create bins by equally spaced location}
}
\value{
A list with bins and values
}
\description{
Using bin data to emperical estimate covariances
}
\details{
Standard empirical covariance estimation.
}
