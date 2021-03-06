% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operator.R
\name{create_operator}
\alias{create_operator}
\title{Create operator components.}
\usage{
create_operator(
  locs,
  name = "matern",
  right.boundary = "neumann",
  left.boundary = "neumann",
  common.grid = FALSE,
  extend = NULL,
  max.dist,
  cutoff = 1e-10,
  n.cores = 1
)
}
\arguments{
\item{locs}{A numeric list of measurement locations.}

\item{name}{A character string for the operator type,
possible options are \code{"matern"}, \code{"exponential"} and \code{"fd2"}.}

\item{right.boundary}{A character string denoting the boundary condition
for the right boundary.}

\item{left.boundary}{A character string denoting the boundary condition
for the left boundary.}

\item{common.grid}{A logical variable for using a common grid for all subjects,
\code{"TRUE"} indicates using a common grid,
\code{"FALSE"} uncommon grids.}

\item{extend}{A numeric vector with two elements specifying the amount of extension
of the grid to the left and right beyondthe measurement locations.}

\item{max.dist}{A numeric value for largest distance between nodes in the discretization}
}
\description{
A function to compute a list of objects for the operator.
}
\details{
This is a supplementary function to be used internally by other functions.
}
\examples{
  \dontrun{
  create_operator(...)
  }

}
\seealso{
\code{\link{estimateLong}}
}
