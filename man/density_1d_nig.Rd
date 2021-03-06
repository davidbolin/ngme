% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CharFunc.R
\name{density_1d_nig}
\alias{density_1d_nig}
\title{Density for NIG stochastic process on an interval}
\usage{
density_1d_nig(interv, n2, param, maternParam)
}
\arguments{
\item{interv}{- start and endpoint to evalute density}

\item{n2}{- number of grid points in the interval}

\item{param}{- (4 x 1)
- (1) delta - location   parameter
- (2) mu    - assymteric parameter
- (3) sigma - scale      parameter
- (4) nu    - shape      parameter}

\item{maternParam}{- (3 x 1)  nu, tau, kappa}
}
\description{
Computes the density for a random field on \eqn{[interv[0], interv[1]]}
}
