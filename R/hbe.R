#' Simulated HBE dataset.
#'
#' @format A dataset containing 76 simulated 2D particle trajectories sampled at 60 Hz.  The four columns of the data are:
#' \describe{
#'   \item{`id`}{Particle ID (1-76).}
#'   \item{`t`}{Time at which position was recorded (in seconds).}
#'   \item{`x`}{Particle position on `x`-axis at time `t` (in microns).}
#'   \item{`y`}{Particle position on `y`-axis at time `t` (in microns).}
#' }
#' @source Hill, D. B., Vasquez, P. A., Mellnik, J., McKinley, S. A., Vose, A., Mu, F., Henderson, A. G., Donaldson, S. H., Alexis, N. E., Boucher, R. C., and Forest, M. G. (2014), "A Biophysical Basis for Mucus Solids Concentration as a Candidate Biomarker for Airways Disease", PLoS ONE(9), e87681.
"hbe"
