% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crouse_msd.R
\name{crouse_msd}
\alias{crouse_msd}
\title{Continuous-spectrum approximation to the Rouse-GLE model.}
\usage{
crouse_msd(t, alpha, tmin, tmax)
}
\arguments{
\item{t}{Vector of time points at which to calculate the MSD.}

\item{alpha}{Subdiffusion coefficient.}

\item{tmin, tmax}{Endpoints of the subdiffusive regime.}
}
\value{
MSD vector.
}
\description{
Continuous-spectrum approximation to the Rouse-GLE model.
}
