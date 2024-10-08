context("Prony-GLE PSD")

## source("subdiff-testfunctions.R")

## # maximum absolute and relative error between two arrays
## max.diff <- function(x1, x2) {
##   c(abs = max(abs(x1-x2)), rel = max(abs(x1-x2)/max(abs(x1), 1e-8)))
## }

ntest <- 20
test_that("PSD FFT == PSD BM+OU", {
  for(ii in 1:ntest) {
    Temp <- runif(1, 290, 300)
    vsigma <- runif(1, 0.01, 0.02)
    K <- sample(10:200, 1)
    alpha <- runif(1, .1, .9)
    tau <- runif(1, 1e-4, 1e-2)
    lambda <- ((1:K)/K)^(1/alpha)/tau
    fseq <- exp(seq(log(1e-6), log(1e6), len = 1000))
    psd1 <- prony_psd_freq(fseq,
                           lambda = lambda, nu = vsigma, Temp = Temp)
    ## rC <- prony_coeff(lambda)
    psd2 <- prony_psd_time(fseq,
                           lambda = lambda, nu = vsigma, Temp = Temp)
    expect_equal(log(psd1), log(psd2))
  }
})
