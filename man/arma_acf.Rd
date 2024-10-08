% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arma_acf.R
\name{arma_acf}
\alias{arma_acf}
\title{Calculate the autocorrelation of an ARMA(p,q) filter applied to a stationary process.}
\usage{
arma_acf(acf, phi = numeric(), rho = numeric(), m)
}
\arguments{
\item{acf}{Vector of length \code{N + m} of autocorrelations of the original stationary process.}

\item{phi}{Vector of \code{p >= 0} coefficients defining the autoregressive part of the filter (see 'Details').}

\item{rho}{Vector of \code{q >= 0} coefficients defining the moving-average part of the filter (see 'Details').}

\item{m}{Order of the MA approximation to the ARMA(p,q) filter.}
}
\value{
A vector of length \code{N} giving the autocorrelation of the filtered process.
}
\description{
Calculate the autocorrelation of an ARMA(p,q) filter applied to a stationary process.
}
\details{
Let \code{X_n} denote the observation of the original process at time \code{n}, and \code{Y_n} denote the corresponding observation of the filtered process.  The AR(p,q) filter model is defined as

\if{html}{\out{<div class="sourceCode">}}\preformatted{Y_n = X_n + sum_\{i=1\}^p phi_i Y_\{n-i\} + sum_\{j=1\}^q rho_j X_\{n-j\}.
}\if{html}{\out{</div>}}

The autocorrelation of \code{Y_n} can be computed exactly from that of \code{X_n} for pure moving-average filters with \code{p = 0}.  The general ARMA(p,q) filter is first approximated by a moving-average process of order \code{m} of the form

\if{html}{\out{<div class="sourceCode">}}\preformatted{Y[n] = X[n] + sum_\{j=1\}^m psi[j] X[n-j],
}\if{html}{\out{</div>}}

where the coefficients \code{psi} are determined using the method of Brockwell & Davis (1991) implemented in \code{\link[stats:ARMAtoMA]{stats::ARMAtoMA()}}.
}
\references{
Brockwell, P.J. and Davis, R.A. "Time Series: Theory and Methods" (1991).  Springer, New York.  \url{https://doi.org/10.1007/978-1-4899-0004-3}.
}
