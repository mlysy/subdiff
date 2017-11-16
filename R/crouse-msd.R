#' Continuous-spectrum approximation to the Rouse-GLE model.
#'
#' @param t Vector of time points at which to calculate the MSD.
#' @param alpha Subdiffusion coefficient.
#' @param tmin,tmax Endpoints of the subdiffusive regime.
#' @return MSD vector.
#' @export
crouse.msd <- function(t, alpha, tmin, tmax) {
  s <- log(t)
  # subdiffusive range
  smin <- log(tmin)
  smax <-log(tmax)
  srng <- smax-smin
  iL <- s < smin # outer left
  iU <- s > smax # outer right
  ind <- !(iL | iU) # inner, i.e., subdiffusion
  # calculate msd
  a <- (1-alpha) * sqrt(pi/8)*srng #a = 1/a;
  msd <- rep(NA, length(s)) #zeros(size(s));
  if(any(ind)) {
    z <- (s[ind]-smin)/srng
    msd[ind] <- exp(s[ind] - a * pnorm(logit(z)))
  }
  msd[iL] <- t[iL]
  msd[iU] <- exp(-a) * t[iU]
  msd
}

logit <- function(x) log(x) - log(1-x)
