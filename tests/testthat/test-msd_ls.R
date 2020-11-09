
context("msd_ls")

test_that(".weighted_ls is equivalent to lm with 'weights'", {
  ntest <- 2
  for(ii in 1:ntest) {
    # simulate data
    alpha <- runif(1, 0, 2)
    D <- rexp(1)
    N <- sample(10:2000, 1)
    ww <- abs(rnorm(N))+1
    xx <- rnorm(N)
    yy <- log(D) + alpha * xx + rexp(1) * rnorm(N)
    # lm method
    theta_lm <- coef(lm(yy ~ xx, weights = ww))
    theta_lm[1] <- exp(theta_lm[1])
    theta_lm <-  theta_lm[2:1]
    names(theta_lm) <- NULL
    # custom method
    theta_w <- .weighted_ls(yy, xx, ww/sum(ww))
    expect_equal(theta_lm, theta_w)
  }
})

test_that("msd_ls pools estimators properly", {
  test_cases <- expand.grid(logw = c(TRUE, FALSE),
                            pooled = c(TRUE, FALSE))
  ntest <- nrow(test_cases)
  for(ii in 1:ntest) {
    # msd conditions
    logw <- test_cases$logw[ii]
    pooled <- test_cases$pooled[ii]
    # simulate data
    alpha <- runif(1, 0, 2)
    D <- rexp(1)
    ntimes <- sample(10:2000, 1)
    npaths <- sample(10:20, 1)
    dT <- runif(1)
    tseq <- (1:ntimes) * dT
    msd <- exp(log(D) + alpha * log(tseq) +
               matrix(rnorm(ntimes*npaths), ntimes, npaths))
    # add NA's
    msd[rbinom(ntimes*npaths, 1, prob = .2) == 1] <- NA
    # msd_ls estimate
    aD_ls <- msd_ls(msd, tseq = tseq, pooled = pooled, logw = logw)
    # lm estimate
    yy <- log(msd)
    xx <- matrix(log(tseq), ntimes, npaths)
    ww <- if(logw) 1/tseq else rep(1, ntimes)
    ww <- matrix(ww, ntimes, npaths)
    if(pooled) {
      aD_lm <- coef(lm(c(yy) ~ c(xx), weights = c(ww)))
      aD_lm <- setNames(c(aD_lm[2], exp(aD_lm[1])), c("alpha", "D"))
    } else {
      aD_lm <- matrix(NA, 2, npaths)
      for(ii in 1:npaths) {
        aD_lm[,ii] <- coef(lm(yy[,ii] ~ xx[,ii], weights = ww[,ii]))
      }
      aD_lm <-  rbind(alpha = aD_lm[2,], D = exp(aD_lm[1,]))
    }
    expect_equal(aD_ls, aD_lm)
  }
})
