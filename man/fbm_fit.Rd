% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fbm_fit.R
\name{fbm_fit}
\alias{fbm_fit}
\title{Fit the fBM model.}
\usage{
fbm_fit(
  Xt,
  dt,
  drift = c("linear", "none", "quadratic"),
  vcov = TRUE,
  ad_only = TRUE
)
}
\arguments{
\item{Xt}{Matrix of trajectory positions, where each row is an observation and each column is a measurement coordinate.  The trajectory is assumed to be sampled at a constant frequency.}

\item{dt}{Interobservation time \eqn{\Delta t} = 1/fps (positive scalar).}

\item{drift}{Character string specifying one of the \emph{preset} drift types "linear", "none", or "quadratic".  Custom drift functions are not supported through this simplified interface.  See \link{csi_model} for details.}

\item{vcov}{Logical; if \code{TRUE}, also estimate the variance matrix.}

\item{ad_only}{Whether to return estimates of (\verb{alpha, log(D)}) only, or the entire parameter vector \code{omega} in the computational basis.}
}
\value{
A vector of estimated parameters on the transformed scale. If \code{vcov == TRUE}, a list with components:
\describe{
\item{coef}{A vector of estimated parameters on transformed scale.}
\item{vcov}{A matrix of estimated covariance of parameters on transformed scale.}
}
If \code{ad_only == TRUE}, instead of the transformed scale parameters, returns an estimate (and possibly the estimated convariance) of \verb{(alpha, D)}.
}
\description{
Fit the fBM model.
}
\examples{
# simulate data from the fbm model
alpha <- .8

dt <- 1/60
N <- 1800
ndim <- 2

Xt <- csi_sim(drift = matrix(0, N-1, ndim),
              acf = fbm_acf(alpha, dt, N-1),
              Sigma = diag(ndim),
              X0 = rep(0, ndim))


# Fit the fbm model
fbm_fit(Xt, dt = dt, drift = "linear")

}
\seealso{
\link{fbm_model}, the class definition for the fBM model.
}
