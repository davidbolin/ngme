% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.functions.R
\name{plot.predict.ngme}
\alias{plot.predict.ngme}
\title{Prediction plots.}
\usage{
\method{plot}{predict.ngme}(object, id = NULL)
}
\arguments{
\item{object}{A fitted object returned by the \code{"predict.ngme"} function.}

\item{id}{A numerical value or character string for ID of the subject
for whom the plot will be generated.}

\item{plot_excursions}{A logical for plotting excursions.}

\item{col_m}{A character value for defining the colour of prediction mean
or median.}

\item{col_c}{A character value for defining the colour of prediction intervals.}

\item{col_p}{A character value for defining the colour of observed data.}

\item{...}{Additional arguments; none used currently.}
}
\description{
Plots the predicted values for a specific subject.
}
\examples{
  \dontrun{
  fit <- ngme(...)
  pred <- predict(fit, ...)
  plot(pred, 1)
  }
}
\seealso{
\code{\link{predict.ngme}}
}
