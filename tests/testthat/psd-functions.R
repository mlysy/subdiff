#--- spectral density of Prony-GLE ---------------------------------------------

#' PSD for the Prony-GLE computed directly in frequency space.
#'
#' @param fseq Frequency vector of length `N`.
#' @param lambda Vector of length `K` giving the inverse decorrelation times of the exponential decay terms.
#' @param nu Scale factor.
#' @param Temp Temperature (Kelvin).
#' @return PSD vector of length `N`.
prony_psd_freq <- function(fseq, lambda, nu = 1, Temp = 298) {
  kB <- 1.3806488e-23 # Boltzmann constant, J/K
  K <- length(lambda)
  N <- length(fseq)
  Sf <- rep(0, N) # force PSD
  G <- rep(0, N) # half-line memory kernel
  for(ii in 1:K){
    Sf <- Sf + lambda[ii] / (lambda[ii]^2 + (2 * pi * fseq)^2)
    G <- G + 1 / (lambda[ii] + 2i * pi * fseq)
  }
  Sf <- 2 * (nu*kB*Temp) * Sf # PSD of laplace distr has factor of 2
  # position scale
  Sx <-  Sf / abs(nu * (2i * pi * fseq) * G)^2
  Sx
}

#' PSD for the Prony-GLE taking the FFT of the time-domain representation.
#'
#' @param fseq Frequency vector of length `N`.
#' @param lambda Vector of length `K` giving the inverse decorrelation times of the exponential decay terms.
#' @param nu Scale factor.
#' @param Temp Temperature (Kelvin).
#' @return PSD vector of length `N`.
prony_psd_time <- function(fseq, lambda, nu = 1, Temp = 298) {
  kB <- 1.3806488e-23 # Boltzmann constant, J/K
  K <- length(lambda)
  N <- length(fseq)
  # coefficients
  rC <- subdiff::prony_coeff(lambda, nu)
  r <- rC$r
  C <- rC$C
  Sx <- C[1]^2 / (2 * pi * fseq)^2 # BM component
  # OU components
  for(ii in 1:(K-1)){
    Sx <- Sx + C[ii+1]^2 /(r[ii]^2 + (2 * pi * fseq)^2)
  }
  ## Sx <- 2 * kB * Temp / nu * Sx
  ## Sx
  kB * Temp * Sx # 2/nu already in C
}


## # psd directly in frequency space
## prony_psd_lambda <- function(fseq, lambda, vsigma, Temp = 298) {
##   kB <- 1.3806488e-23 # Boltzmann constant, J/K
##   K <- length(lambda)
##   N <- length(fseq)
##   Sf <- rep(0, N) # force PSD
##   G <- rep(0, N) # half-line memory kernel
##   for(ii in 1:K){
##     Sf <- Sf + lambda[ii] / (lambda[ii]^2 + (2 * pi * fseq)^2)
##     G <- G + 1 / (lambda[ii] + 2i * pi * fseq)
##   }
##   Sf <- 2 * (vsigma*kB*Temp) * Sf # PSD of laplace distr has factor of 2
##   # position scale
##   Sx <-  Sf / abs(vsigma * (2i * pi * fseq) * G)^2
##   Sx
## }

## # psd from r/C time-domain representation
## prony_psd_rC <- function(fseq, lambda, vsigma, Temp) {
##   kB <- 1.3806488e-23 # Boltzmann constant, J/K
##   K <- length(lambda)
##   N <- length(fseq)
##   # coefficients
##   rC <- subdiff::prony_coeff(lambda)
##   r <- rC$r
##   C <- rC$C
##   Sx <- C[1]^2 / (2 * pi * fseq)^2 # BM component
##   # OU components
##   for(ii in 1:(K-1)){
##     Sx <- Sx + C[ii+1]^2 /(r[ii]^2 + (2 * pi * fseq)^2)
##   }
##   Sx <- 2 * kB * Temp / vsigma * Sx
##   Sx
## }
