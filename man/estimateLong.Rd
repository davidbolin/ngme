% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateGH.R
\name{estimateLong}
\alias{estimateLong}
\title{Estimate parameters.}
\usage{
estimateLong(
  Y,
  locs = NULL,
  mixedEffect_list,
  measurment_list,
  processes_list,
  operator_list,
  step0 = 0.3,
  alpha = 0.3,
  learning_rate = 0,
  nBurnin_learningrate = NULL,
  nBurnin_base = 0,
  pSubsample = 1,
  polyak_rate = -1,
  subsample.type = 1,
  iter_start = 0,
  nPar_burnin = 0,
  pSubsample2 = 0.3,
  nIter = 10,
  nSim = 1,
  nBurnin = 10,
  silent = FALSE,
  seed = NULL,
  standardize.mixedEffects = FALSE,
  estimate_fisher = FALSE
)
}
\arguments{
\item{Y}{A numeric list that contains outcome values.}

\item{locs}{A numeric list that contains the timings at which the outcomes
are collected.}

\item{mixedEffect_list}{A list of inputs for random effects.
\itemize{
\item \code{noise} The distribution of the mixed effects.
\item \code{B_random} A list that contains the random effect
   covariates (needs to be matrix, can be NULL).
\item \code{B_fixed} A list that contains the fixed effect
   covariates (needs to be matrix, can be NULL).
\item \code{beta_random} Initial values for the parameters of the
   random effects (mean parameter) (if not specified set to zero).
\item \code{beta_fixed} Initial  values for the parameters of the
   fixed effects (if not specified set to zero).
\item \code{Sigma} Initial values for the parameters of the
   variance-covariance matrix of the random effects
   (if not specified set to I ).
\item \code{nu} Shape parameter for noise (NIG only)
\item \code{mu} Shift parameter for noise (NIG only)
\item \code{U} A list of inital values of the random effects.
\item \code{V} A list of inital values of the variance effects.
}}

\item{measurment_list}{A list of inputs for measurement error.
\itemize{
\item \code{sigma} Measurement noise variance parameter.
\item \code{nu}    Shape parameter for noise (NIG only).
\item \code{Vs}    A list of inital values for the noise of the measurement.
}}

\item{processes_list}{A list of inputs for the process.
\itemize{
\item \code{noise} Distribution of the process.
\item \code{nu}    Shape parameter (for NIG or GAL).
\item \code{mu}    Asymmetry parameter (for NIG or GAL).
}}

\item{step0}{A numeric value for stepsize for the optimizer; step0 / i^alpha.}

\item{alpha}{A numeric value for stepsize for the optimizer; step0 / i^alpha.}

\item{learning_rate}{A numeric value for the parameter of stochastic gradient.}

\item{nBurnin_learningrate}{A numeric value until which the learning will
not be started.}

\item{nBurnin_base}{A numerical value for burn-in simulations that are performed
for a subject that is sampled for the first time in the estimation method.}

\item{pSubsample}{A numeric value for the portion of data to be used in each
gradient iteration.}

\item{polyak_rate}{A numeric value for moving average of parameters;
-1: inactive, 0: pure mean.}

\item{subsample.type}{A numeric value for the type of subsampling;
1: uniform sampling,
2: sample size weighted,
3: weighted sampling by gradient size,
4: grouped sub-sampler.}

\item{nPar_burnin}{A numeric value; "M-step" updates will be used until this
iteration.}

\item{pSubsample2}{A numeric value for the portion of the data
to be used in each gradient subsampling weighted by gradient.}

\item{nIter}{A numeric value for the number of iteration that will be
used by the stochastic gradient.}

\item{nSim}{A numeric value for the number of samples of the Gibbs sampler
to estimate the gradient.}

\item{silent}{A logical value for printing the details of the iterations;
\code{"TRUE"} indicates do not print, \code{"FALSE"} indicates print.}

\item{seed}{A numerical value for starting the Gibbs samplers from fixed seed.}

\item{standardize.mixedEffects}{A logical variable for standardising the covariates;
\code{"FALSE"} indicates no standardisation, \code{"TRUE"} standardisation.}

\item{estimate_fisher}{A logical variable for whether Fisher-Information matrix
to be obtained; \code{"FALSE"} indicates do not obtain, \code{"TRUE"} obtain.}
}
\value{
A list of output.
}
\description{
A function that estimates parameters by
   calling the \code{"estimateLong_cpp()"} function.
}
\details{
This function calls \code{"estimateLong_cpp()"} internally.
   It is wrapped by \code{"ngme"}, and is not advised to
   be used alone.
}
\examples{
  \dontrun{
  data(srft_data)
  fit <- estimateLong(...)
  }
}
\seealso{
\code{\link{ngme}}
}
