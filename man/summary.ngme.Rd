% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngme.utils.R
\name{summary.ngme}
\alias{summary.ngme}
\title{Summary function for \code{"ngme"} objects.}
\usage{
\method{summary}{ngme}(object, ...)
}
\arguments{
\item{object}{A fitted object by calling \code{"ngme"}.}

\item{...}{Additional arguments; none used currently.}
}
\description{
A function to summarise the output contained in 
   \code{"ngme"} objects.
}
\examples{
  \dontrun{
  fit <- ngme(...)
  summary(fit)
  }
}
\seealso{
\code{\link{ngme}}, \code{\link{print.summary.ngme}},
  \code{\link{print}}, \code{\link{summary}}
}
