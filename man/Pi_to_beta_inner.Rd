% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pi_to_beta_inner.R
\name{Pi_to_beta_inner}
\alias{Pi_to_beta_inner}
\title{Maps the joint distribution, Pi, of X and Y* to a scalar, beta}
\usage{
Pi_to_beta_inner(Pi, X_vals, Y_vals, W_weights)
}
\arguments{
\item{Pi}{A numeric vector or list of numeric vectors containing the elements of Pi.}

\item{X_vals}{A numeric vector or a list of numeric vectors representing the scalar values associated with X.}

\item{Y_vals}{A numeric vector or a list of numeric vectors representing the scalar values associated with Y.}

\item{W_weights}{A numeric vector representing the sample size of each control cell.}
}
\value{
A scalar equal to beta.
}
\description{
Maps the joint distribution, Pi, of X and Y* to a scalar, beta
}
\keyword{internal}
