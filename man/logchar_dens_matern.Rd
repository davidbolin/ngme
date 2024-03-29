% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CharFunc.R
\name{logchar_dens_matern}
\alias{logchar_dens_matern}
\title{log-characterisct function of integral}
\usage{
logchar_dens_matern(t, logphi_W, maternParam, d, typeMatern = 1)
}
\description{
log charateristic function of \eqn{\int f(x)dW_x} assuming that f is istorpic,
where W_x brownian sheet.

t           - (vector) location of chara eval
logphi_W    - (function)   evaluates log of char function of \eqn{f(x)*W(dx)}
maternParam - (3 x 1)  nu, tau, kappa
d           - (1 x 1)  dimension
}
