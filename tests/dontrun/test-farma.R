#--- check speed of the farma_acf calculations ---------------------------------

require(subdiff)
require(fftw)

N <- 2000
m <- 30
nreps <- 100

# fma
system.time({
  replicate(nreps, {
    farma_acf(alpha = runif(1), phi = 0, rho = runif(3, -.5, .5),
              dT = runif(1), N = N)
  })
})

# far
system.time({
  replicate(nreps, {
    farma_acf(alpha = runif(1), phi = runif(1, -.5, .5), rho = 0,
              dT = runif(1), N = N, m = m)
  })
})

# farma
system.time({
  replicate(nreps, {
    farma_acf(alpha = runif(1), phi = runif(1, -.5, .5), rho = runif(1, -.5, .5),
              dT = runif(1), N = N, m = m)
  })
})


system.time({
  plan <- planFFT(N+m)
  replicate(nreps, {
    FFT(rnorm(N+m), plan = plan)
    FFT(rnorm(N+m), plan = plan)
    FFT(rnorm(N+m), plan = plan)
    FFT(rnorm(N+m), plan = plan)
    FFT(rnorm(N+m), plan = plan)
    FFT(rnorm(N+m), plan = plan)
    FFT(rnorm(N+m), plan = plan)
    FFT(rnorm(N+m), plan = plan)
  })
})

system.time({
  replicate(nreps, {
    FFT(rnorm(N+m))
    FFT(rnorm(N+m))
    FFT(rnorm(N+m))
    FFT(rnorm(N+m))
    FFT(rnorm(N+m))
    FFT(rnorm(N+m))
    FFT(rnorm(N+m))
    FFT(rnorm(N+m))
  })
})
