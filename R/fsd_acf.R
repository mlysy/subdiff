#' Calculate the autocorrelation of the fSD model.
#'
#' Compute the autocorrelation of the Savin & Doyle (2005) localization error model with fBM increments (see 'Details').
#'
#' @param alpha Subdiffusion exponent of the underlying fBM process. A scalar between 0 and 2.
#' @param tau The ratio between camera shutter open time and interobservation time `dt` (see 'Details'). A scalar between 0 and 1.
#' @param sigma2 The magnitude of the static error (see 'Details'). A positive scalar.
#' @template args-dt
#' @template args-N
#'
#' @template ret-acf
#'
#' @details Let `X_t` denote the position of an fBM process at time `t`.  The Savin-Doyle localization error model describes the observed position `Y_n` at time `t = n * dt` as
#' ```
#' Y_n = sigma * eps_n + 1/tau * int_0^tau X_{n*dt + s} ds,
#' ```
#' where `eps_n ~iid N(0,1)` is a Gaussian white noise process.
#'
#' This function returns the autocorrelation of the stationary process `dY_n = Y_{n+1} - Y_n`.
#'
#' @references Savin, T., and Doyle, P.S. "Static and dynamic errors in particle tracking microrheology." Biophysical Journal 88.1 (2005): 623-638.
#'
#' @export
fsd_acf <- function(alpha, tau, sigma2, dt, N) {
  if(tau == 0) {
    acf1 <- fbm_acf(alpha, dt, N)
  } else {
    vec <- (0:N) / tau
    vec <- (vec + 1)^(alpha+2) + abs(vec - 1)^(alpha+2) - 2 * vec^(alpha+2)
    if(N == 1) {
      acf1 <- 2 * (vec[2] - vec[1])
    } else {
      acf1 <- vec[1:N + 1] + c(vec[2], vec[1:(N - 1)]) - 2 * vec[1:N]
    }
    acf1 <- .5 * (tau * dt)^alpha / ((alpha+1) * (alpha+2)) * acf1
    ## vec <- .g_func(alpha, tau, dt, N)
    ## if(N == 1) {
    ##   acf1 <- 2 * vec[2] - 2 * vec[1]
    ## } else {
    ##   acf1 <- vec[1:N + 1] + c(vec[2], vec[1:(N - 1)]) - 2 * vec[1:N]
    ## }
  }
  if(N == 1) {
    acf1 <- acf1 + 2 * sigma2
  } else {
    acf1[1:2] <- acf1[1:2] + sigma2*c(2,-1)
  }
  acf1
}

## .g_func <- function(alpha, tau, dt, N) {
##   vec <- (0:N) / tau
##   alpha2 <- alpha + 2
##   ans <- (vec + 1)^alpha2 + abs(vec - 1)^alpha2 - 2 * vec^alpha2
##   ans * (tau * dt)^alpha / (alpha+1) / alpha2 / 2
## }

#--- quick test ----------------------------------------------------------------

## alpha <- .3
## dt <- .11
## N <- 5
## tau <- .46

## fsd_acf(alpha = alpha, tau = tau, sigma2 = 0, dt = dt, N = N) - fbm_acf(alpha = alpha + 2, dt = 1/tau, N = N)

## alpha <- runif(1, 0, 2)
## dt <- runif(1)
## N <- sample(5, 1)
## tau <- runif(1)
## sigma2 <- runif(1)
## fsd_acf(alpha, tau, sigma2, dt, N) - fsd_acf2(alpha, tau, sigma2, dt, N)
