# spectral density of Prony-GLE

# psd directly in frequency space
prony.psd.lambda <- function(fseq, lambda, vsigma, Temp = 298) {
  kB <- 1.3806488e-23 # Boltzmann constant, J/K
  K <- length(lambda)
  N <- length(fseq)
  Sf <- rep(0, N) # force PSD
  G <- rep(0, N) # half-line memory kernel
  for(ii in 1:K){
    Sf <- Sf + lambda[ii] / (lambda[ii]^2 + (2 * pi * fseq)^2)
    G <- G + 1 / (lambda[ii] + 2i * pi * fseq)
  }
  Sf <- 2 * (vsigma*kB*Temp) * Sf # PSD of laplace distr has factor of 2
  # position scale
  Sx <-  Sf / abs(vsigma * (2i * pi * fseq) * G)^2
  Sx
}

# psd from r/C time-domain representation
prony.psd.rC <- function(fseq, lambda, vsigma, Temp) {
  kB <- 1.3806488e-23 # Boltzmann constant, J/K
  K <- length(lambda)
  N <- length(fseq)
  # coefficients
  rC <- prony.coeff(lambda)
  r <- rC$r
  C <- rC$C
  Sx <- C[1]^2 / (2 * pi * fseq)^2 # BM component
  # OU components
  for(ii in 1:(K-1)){
    Sx <- Sx + C[ii+1]^2 /(r[ii]^2 + (2 * pi * fseq)^2)
  }
  Sx <- 2 * kB * Temp / vsigma * Sx
  Sx
}
