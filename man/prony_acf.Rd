% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prony_acf.R
\name{prony_acf}
\alias{prony_acf}
\title{Autocorrelation of the Prony-GLE increment process.}
\usage{
prony_acf(lambda, nu = 1, dt, N, ...)
}
\arguments{
\item{lambda}{Vector of length \code{K} giving the inverse decorrelation times of the exponential decay terms.}

\item{nu}{Autocorrelation scaling factor.}

\item{dt}{Interobservation time \eqn{\Delta t} = 1/fps (positive scalar).}

\item{N}{Number of observations (positive integer).}

\item{...}{Additional arguments to \code{\link[=prony_coeff]{prony_coeff()}}.}
}
\value{
A vector of \code{N} autocorrelation values.
}
\description{
Autocorrelation of the Prony-GLE increment process.
}
\details{
The Prony-GLE model satisfies the integro-differential equation
\deqn{
F_t - \int_{-\infty}^t \kappa(t-s) \dot X_s ds = 0, \qquad ACF_F(t) = k_B T \kappa(t),
}
where \eqn{T} is temperature, \eqn{k_B} is Boltzmann's constant, and the memory kernel is a sum of exponentials:
\deqn{
\kappa(t) = \nu sum_{k=1}^K exp(-\lambda_k t).
}
The solution process is of the form
\deqn{
X_t = C_0 B_t + \sum_{i=1}^{K-1} C_i W_{it},
}
where \eqn{B_t} is Brownian motion and \eqn{d W_{it} = -\rho_i W_{it} dt + d B_{it}} are Ornstein-Uhlenbeck processes all independent of each other and of \eqn{B_t}.
}
\note{
Since temperature is not provided, the result is only proportional to the desired ACF, such that

\if{html}{\out{<div class="sourceCode">}}\preformatted{ACF_dX(t) = k_B T * prony_acf(lambda, nu, N, dt)
}\if{html}{\out{</div>}}
}
