% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rouse_sub.R
\name{rouse_sub}
\alias{rouse_sub}
\title{Subdiffusive timescale of Rouse-GLE model.}
\usage{
rouse_sub(alpha, tau, K, est_alpha = TRUE, ...)
}
\arguments{
\item{alpha}{Rouse-GLE subdiffusion coefficient.}

\item{tau}{Shortest timescale of force memory.}

\item{K}{Number of modes in relaxation spectrum.}

\item{est_alpha}{If \code{FALSE}, use the Rouse-GLE's \code{alpha} parameter, otherwise estimate \code{alpha} by least-squares.}

\item{...}{Additional arguments to pass to \code{prony_coeff}.}
}
\value{
Vector with named elements \code{tmin}, \code{tmax}, \code{alpha}, and \code{D}.
}
\description{
Calculates the range of the subdiffusive timescale along with the effective subdiffusion coefficient and diffusivity constant on this range.
}
\details{
The Rouse model is a zero-mass, zero-potential GLE with force memory
kernel
\deqn{
acf_F(t) = 1/K sum( exp(- t * (1:K/K)^rho / tau) ).
}
As \eqn{K} gets large, it exhibits has transient subdiffusion:
\deqn{
MSD(t) ~ D * t^alpha, t \in tScale
         t,           t \notin tScale,
}
where \eqn{alpha = 1/rho} is the subdiffusion exponent and \eqn{D} is the diffusivity
constant.
}
