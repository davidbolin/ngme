% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CharFunc.R
\name{maternkernelMulti}
\alias{maternkernelMulti}
\title{kernel function producing the Matern covariance}
\usage{
maternkernelMulti(x, alpha, tau, kappa, d)
}
\arguments{
\item{x}{- (n x 1) distance between points}

\item{tau}{- (>0)    precision parameter}

\item{kappa}{- (>0)    range-scale parameter}

\item{d}{- (int>0) dimension of x}

\item{nu}{- (>0)    smoothnes}
}
\value{
f(x) - (n x 1) kernel evaluted at x
}
\description{
kernel function producing the Matern covariance
}
