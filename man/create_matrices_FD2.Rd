% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operator.R
\name{create_matrices_FD2}
\alias{create_matrices_FD2}
\title{Create matrices for Finite difference.}
\usage{
create_matrices_FD2(
  locs,
  common.grid = FALSE,
  extend = NULL,
  max.dist,
  cutoff = 1e-10,
  n.cores = 1,
  Y = NULL,
  loc.Y = NULL,
  max.dY = -1,
  nJump = 3
)
}
\arguments{
\item{locs}{A numeric list for the meansurement locations.}

\item{common.grid}{A logical variable for using a common grid for all subjects,
\code{"TRUE"} indicates using a common grid,
\code{"FALSE"} uncommon grids.}

\item{extend}{A numeric vector with two elements specifying the amount of extension
of the grid to the left and right beyondthe measurement locations.}

\item{max.dist}{A numeric value for largest distance between nodes in the discretization}
}
\value{
Returns matrices.
}
\description{
A function to create matrices for Finite difference operator, one sided.
}
\details{
This is a supplementary function to be used internally by other functions.
}
\examples{
  \dontrun{
  create_matrices_FD2(...)
  }

}
\seealso{
\code{\link{create_operator}}
}
