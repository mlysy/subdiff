# functions for msd calculations

source(system.file("proj", "plot-functions.R", package = "subdiff"))

#-------------------------------------------------------------------------------

# maximum relative error to ground truth
# msd0: vector
# msd: vector or matrix
max_re <- function(msd0, msd) max(exp(abs(log(msd0/msd))))

# mu_hat is a vector of point estimates
# mu_boot is the bootstrap sample (each column is one draw)
# calculates the basic bootstrap
boot_conf <- function(mu_hat, mu_boot, conf = .95) {
  conf <- (1-conf)/2
  ci <- apply(2*mu_hat - mu_boot, 1, quantile, probs = c(conf, 1-conf))
  rownames(ci) <- c("L", "U")
  cbind(est = mu_hat, t(ci))
}

# msd est + ci from empirical MSDs.
emp_conf <- function(msd, nboot = 1000) {
  npaths <- ncol(msd)
  mu_hat <- rowMeans(msd, na.rm = TRUE)
  # simple bootstrap confidence interval
  mu_boot <- replicate(nboot, {
    ind <- sample(npaths, npaths, replace = TRUE)
    rowMeans(msd[,ind,drop=FALSE], na.rm = TRUE)
  })
  boot_conf(mu_hat, mu_boot)
  ## ci <- apply(mu.boot, 1, quantile, probs = c(.025, .975))
  ## rownames(ci) <- c("L", "U")
  ## cbind(est = mu, t(ci))
}

# msd est + ve + ci by ls method
# fit_mean: apply ls estimator to overall mean, instead of individual trajectories
ls_conf <- function(msd, tfit, tseq, logw = FALSE, fit_mean = FALSE, nboot = 1000) {
  npaths <- ncol(msd)
  if(!fit_mean) {
    # point estimate
    aD_hat <- msd_ls(msd = msd, tseq = tfit, pooled = TRUE, logw = logw)
    # bootstrap
    aD_boot <- replicate(nboot, {
      msd_ls(msd = msd[,sample(npaths, replace = TRUE)],
             tseq = tfit, pooled = TRUE, logw = logw)
    })
  } else {
    # first bootstrap the mean estimate
    mu_hat <- rowMeans(msd, na.rm = TRUE)
    mu_boot <- replicate(nboot, {
      ind <- sample(npaths, npaths, replace = TRUE)
      rowMeans(msd[,ind,drop=FALSE], na.rm = TRUE)
    })
    # now extract aD estimates
    aD_hat <- msd_ls(mu_hat, tseq = tfit, pooled = TRUE, logw = logw)
    aD_boot <- msd_ls(mu_boot, tseq = tfit, pooled = FALSE, logw = logw)
  }
  # confidence interval
  aD_conf(aD_hat, t(aD_boot), tseq)
  ## list(aD_hat = aD_hat, aD_boot = t(aD_boot))
}

# msd est + ve + ci by fbm method
fbm_conf <- function(alpha, logD, tseq, nboot = 1000) {
  aD <- cbind(alpha, logD)
  npaths <- nrow(aD)
  # point estimate
  aD_hat <- colMeans(aD)
  aD_hat[2] <- exp(aD_hat[2])
  names(aD_hat)[2] <- "D"
  # bootstrap
  aD_boot <- replicate(nboot, {
    colMeans(aD[sample(npaths, npaths, replace = TRUE),])
  })
  aD_boot <- cbind(alpha = aD_boot[1,], D = exp(aD_boot[2,]))
  aD_conf(aD_hat, aD_boot, tseq)
}

## fbm_conf_old <- function(alpha, D, tseq, nboot = 1000) {
##   npaths <- length(alpha)
##   logD <- log(D)
##   aD_hat <- c(alpha = mean(alpha), D = exp(mean(logD)))
##   # bootstrap
##   aD_boot <- replicate(nboot, {
##     bind <- sample(npaths, npaths, replace = TRUE)
##     c(alpha = mean(alpha[bind]), D = exp(mean(logD[bind])))
##   })
##   # confidence interval
##   aD_conf(aD_hat, t(aD_boot), tseq)
## }

# msd est + ci based on point and bootstap aD = (alpha, D) pairs
# also outputs coef = (alpha, logD) and vcov = variance estimate of these
aD_conf <- function(aD_hat, aD_boot, tseq) {
  aD_hat[2] <- log(aD_hat[2])
  aD_boot[,2] <- log(aD_boot[,2])
  mu_boot <- apply(aD_boot, 1, function(theta) {
    theta[2] + theta[1] * log(tseq)
  })
  mu_hat <- aD_hat[2] + aD_hat[1] * log(tseq)
  list(coef = aD_hat, vcov = var(aD_boot),
       msd = exp(boot_conf(mu_hat, mu_boot)))
}

#--- msd_plots -----------------------------------------------------------------

# extract path with given id
get_path <- function(pathid, paths, names = c("X", "Y")) {
  as.matrix(paths[paths$pathID == pathid, names])
}

msd_ci <- function(msd, npaths = ncol(msd), nboot = 1000) {
  mu <- rowMeans(msd, na.rm = TRUE)
  # delta method for ci on lambda = log(mu)
  ## n <- ncol(x)
  mu.boot <- replicate(nboot, {
    ind <- sample(ncol(msd), npaths, replace = TRUE)
    log(rowMeans(msd[,ind,drop=FALSE], na.rm = TRUE))
  })
  sig <- apply(mu.boot, 1, sd)
  ci <- exp(log(mu) + cbind(L = -2*sig, U = 2*sig))
  cbind(est = mu, ci)
}

get_theta <- function(pathid, stats) {
  theta <- pstats$coef
  setNames(nm = c("alpha", "logD"),
           object = c(itrans_alpha(theta["gamma"]), theta["lambda1"]))
}

# distance between empirical and theoretical msd
# each column is a realization.
# msd_theo can be a vector or matrix, in which case average is first taken
# at each time point.
msd_dist <- function(msd, alpha, logD, dT, msd_theo,
                     logw = TRUE, bytime = FALSE) {
  nlag <- nrow(msd)
  ww <- if(logw) 1/(1:nlag) else rep(1, nlag)
  ww <- ww/sum(ww)
  if(missing(msd_theo)) {
    msd_theo <- logD + alpha * log((1:nlag-1)*dT)
  } else {
    msd_theo <- log(rowSums(msd_theo))
  }
  mdist <- ww * (log(msd) - msd_theo)^2
  if(bytime) rowSums(mdist) else sum(mdist)
}

# msd per-trajectory + various estimates of overall msd
msd_plot <- function(tseq, ci_list, msd, ci_clrs, msd_clrs, xlim, ord) {
  # plot msd
  if(missing(xlim)) xlim <- range(tseq)
  tind <- (tseq <= xlim[2]) & (tseq >= xlim[1]) # corresponding time points
  ylim <- range(do.call(cbind, c(list(msd = msd), ci_list))[tind,],
                na.rm = TRUE)
  .multi_plot(tseq, msd, log = "xy", xlim, ylim = ylim, col = msd_clrs)
  # plot CI's
  nci <- length(ci_list)
  if(missing(ord)) ord <- 1:nci
  for(ii in ord) {
    .multi_plot(tseq, ci_list[[ii]], lty = c(1,2,2), lwd = c(2,1,1),
               col = ci_clrs[ii], add = TRUE)
  }
  invisible(NULL)
}


