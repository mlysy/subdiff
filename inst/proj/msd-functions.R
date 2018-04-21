# functions for msd calculations

source(system.file("proj", "plot-functions.R", package = "subdiff"))

#-------------------------------------------------------------------------------

# maximum relative error to ground truth
# msd0: vector
# msd: vector or matrix
max_re <- function(msd0, msd) {
  max(exp(abs(log(msd0/msd))))
  ## max(abs(1-msd/msd0))
}

# mu_hat is a vector of point estimates
# mu_boot is the bootstrap sample (each column is one draw)
# calculates the basic bootstrap
boot_ci <- function(mu_hat, mu_boot, conf = .95) {
  conf <- (1-conf)/2
  ci <- apply(2*mu_hat - mu_boot, 1, quantile, probs = c(conf, 1-conf))
  rownames(ci) <- c("L", "U")
  cbind(est = mu_hat, t(ci))
}

# msd est + ci from empirical MSDs.
emp_msdCI <- function(msd, nboot = 1000) {
  npaths <- ncol(msd)
  mu_hat <- rowMeans(msd, na.rm = TRUE)
  # simple bootstrap confidence interval
  mu_boot <- replicate(nboot, {
    ind <- sample(npaths, npaths, replace = TRUE)
    rowMeans(msd[,ind,drop=FALSE], na.rm = TRUE)
  })
  boot_ci(mu_hat, mu_boot)
  ## ci <- apply(mu.boot, 1, quantile, probs = c(.025, .975))
  ## rownames(ci) <- c("L", "U")
  ## cbind(est = mu, t(ci))
}

# msd est + ci by ls method
ls_msdCI <- function(msd, tfit, tseq, logw = FALSE, nboot = 1000) {
  npaths <- ncol(msd)
  # point estimate
  aD_hat <- msd_ls(msd = msd, tseq = tfit, pooled = TRUE, logw = logw)
  # bootstrap
  aD_boot <- replicate(nboot, {
    msd_ls(msd = msd[,sample(npaths, replace = TRUE)],
           tseq = tfit, pooled = TRUE, logw = logw)
  })
  # confidence interval
  aD_msdCI(aD_hat, t(aD_boot), tseq)
  ## list(aD_hat = aD_hat, aD_boot = t(aD_boot))
}

fbm_msdCI <- function(alpha, D, tseq, nboot = 1000) {
  npaths <- length(alpha)
  logD <- log(D)
  aD_hat <- c(alpha = mean(alpha), D = exp(mean(logD)))
  # bootstrap
  aD_boot <- replicate(nboot, {
    bind <- sample(npaths, npaths, replace = TRUE)
    c(alpha = mean(alpha[bind]), D = exp(mean(logD[bind])))
  })
  # confidence interval
  aD_msdCI(aD_hat, t(aD_boot), tseq)
}

# msd est + ci based on point and bootstap aD = (alpha, D) pairs
aD_msdCI <- function(aD_hat, aD_boot, tseq) {
  aD_boot[,2] <- log(aD_boot[,2])
  mu_boot <- apply(aD_boot, 1, function(theta) {
    theta[2] + theta[1] * log(tseq)
  })
  mu_hat <-  log(aD_hat[2]) + aD_hat[1] * log(tseq)
  exp(boot_ci(mu_hat, mu_boot))
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
msd_plot <- function(tseq, ci_list, msd, clrs, xlim, ord) {
  # plot msd
  if(missing(xlim)) xlim <- range(tseq)
  tind <- (tseq <= xlim[2]) & (tseq >= xlim[1]) # corresponding time points
  ylim <- range(do.call(cbind, c(list(msd = msd), ci_list))[tind,],
                na.rm = TRUE)
  .multi_plot(tseq, msd, log = "xy", xlim, ylim = ylim, col = clrs[1])
  # plot CI's
  nci <- length(ci_list)
  if(missing(ord)) ord <- 1:nci
  for(ii in ord) {
    .multi_plot(tseq, ci_list[[ii]], lty = c(1,2,2), lwd = c(2,1,1),
               col = clrs[ii+1], add = TRUE)
  }
  invisible(NULL)
}


