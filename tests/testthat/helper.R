#--- subdiff testing helper functions ------------------------------------------

## # set the seed if we're on cran
## if(!identical(Sys.getenv("NOT_CRAN"), "true")) {
##   set.seed(2022)
## }
set.seed(2022) # so randomized tests never fail

# dimension of covariance matrix in computational basis
getq <- function(ndims) {
  ## if(ndims == 1) 1 else 3
  ndims * (ndims+1) / 2
}

# max of min of abs and rel error
max_xdiff <- function(x) {
  xdiff <- abs(diff(x))
  max(pmin(xdiff[,1], xdiff[,2]))
}

# log(1+exp(x))
log1pexp <- function(x) {
  n <- length(x)
  y <- rep(NA, n)
  ind <- x <= -37
  if(any(ind)) y[ind] <- exp(x[ind])
  ind <- x > -37 & x <= 18
  if(any(ind)) y[ind] <- log1p(exp(x[ind]))
  ind <- x > 18 & x <= 33.3
  if(any(ind)) y[ind] <- x[ind] + exp(-x[ind])
  ind <- x > 33.3
  if(any(ind)) y[ind] <- x[ind]
  y
}


## # inverse logit function
## ilogit <- function(x, min = 0, max = 1) {
##   1/(1+exp(-x)) * (max - min) + min
## }

## logit <- function(x, min = 0, max = 1) {
##   x <- (x - min) / (max - min)
##   log(x) - log(1-x)
## }

farma_acf2 <- function(alpha, phi, rho, dt, N, m = 30) {
  # ACF of fbm process (internal)
  .fbm_acf <- function(alpha, dt, N) {
    if(N == 1) {
      acf <- dt^alpha
    } else {
      acf <- (dt*(0:N))^alpha
      acf <- .5 * (acf[1:N+1] + c(acf[2], acf[1:(N-1)]) - 2*acf[1:N])
    }
    acf
  }
  # ACF of unparametrized moving-average model with fBM noises.
  .fma_acf <- function(alpha, rho, dt, N) {
    nlag <- length(rho)
    acf1 <- .fbm_acf(alpha, dt, N+nlag)
    if(nlag == 1) {
      acf2 <- rho^2 * acf1[1:N]
    } else if(nlag == 2) {
      acf2 <- sum(rho^2)*acf1[1:N] + (rho[1]*rho[2])*(acf1[c(2,2:N-1)]+acf1[1:N+1])
    } else if(nlag == 3) {
      acf2 <- sum(rho^2)*acf1[1:N] +
        (rho[1]*rho[2]+rho[2]*rho[3])*(acf1[c(2,2:N-1)]+acf1[1:N+1]) +
        (rho[1]*rho[3])*(acf1[c(3,2,3:N-2)]+acf1[1:N+2])
    } else {
      a <- c(rho[1], rep(0, N), rev(rho[-1]))
      b <- c(acf1, rev(acf1[2:nlag]))
      c <- c(rho, rep(0, N+nlag-1))
      k1 <- Re(fftw::IFFT(fftw::FFT(b)*fftw::FFT(c)))[1:(N+nlag)]
      acf2 <- Re(fftw::IFFT(fftw::FFT(a)*fftw::FFT(k1)))[1:N]
    }
    acf2
  }
  # ACF of un-parametrizd far(1)
  .ar1_acf <- function(acf1, phi, N, m) {
    a <- c(1, rep(0, N), phi^(m:1))
    b <- c(acf1, rev(acf1[1:m+1]))
    c <- c(phi^(0:m), rep(0, N+m))
    k1 <- Re(fftw::IFFT(fftw::FFT(b)*fftw::FFT(c)))[1:(N+m+1)]
    acf2 <- Re(fftw::IFFT(fftw::FFT(a)*fftw::FFT(k1)))[1:N]
    acf2
  }
  # now farma
  if(!phi) {
    acf2 <- .fma_acf(alpha, c(1-sum(rho), rho), dt, N)
  } else {
    acf1 <- .fma_acf(alpha, c(1-sum(phi)-sum(rho), rho), dt, N+m+1)
    acf2 <- .ar1_acf(acf1, phi, N, m)
  }
  acf2
}

ls_var_r <- function(alpha, lags, N) {
  tprod <- 1/sqrt(lags %o% lags)
  t1ov2 <- sqrt(lags %o% (1/lags))
  t2ov1 <- 1/t1ov2
  V <- 0
  for(ii in (-N+1):(N-1)) {
    Vt <- abs(ii * tprod + t1ov2)^alpha
    Vt <- Vt - abs(ii * tprod + t1ov2 - t2ov1)^alpha
    Vt <- Vt - abs(ii * tprod)^alpha
    Vt <- Vt + abs(ii * tprod - t2ov1)^alpha
    V <- V + (1-abs(ii)/N) * Vt^2
  }
  .5 * V/N
}

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
