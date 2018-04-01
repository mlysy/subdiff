#--- max-jump tests ------------------------------------------------------------

# ratio of maximum to median increment
rmax <- function(dX) {
  mu <- apply(dX, 2, median)
  D <- sqrt(colSums((t(dX) - mu)^2))
  max(D)/median(D)
}

# z-score of max increment instead
zmax <- function(dX, log = FALSE) {
  ## mu <- apply(dX, 2, median)
  mu <- colMeans(dX)
  D <- sqrt(colSums((t(dX) - mu)^2)) # increment sizes
  if(log) D <- log(D)
  (max(D) - mean(D))/sd(D)
}

#--- fit paths and gof tests as well -------------------------------------------

# eventually this will just be an option for the fitting itself

resid_gof <- function(Z) {
  nm <- c("Z1", "Z2", "Zc")
  # anderson-darling test
  ad <- setNames(c(ad_pval(Z[,1]), ad_pval(Z[,2]), ad_pval(c(Z))), nm = nm)
  # modified shapiro-wilks
  sw <- setNames(c(sw_pval(Z[,1]), sw_pval(Z[,2]), sw_pval(c(Z))), nm = nm)
  # berkowitz correlation test
  bc <- c(bc_pval(Z[,1]), bc_pval(Z[,2]))
  bc <- setNames(c(bc, pbeta(min(bc), 1, 2)), nm = nm)
  list(ad = ad, sw = sw, bc = bc)
}

#--- experimentwise pvalue plots -----------------------------------------------

pval_plot <- function(dstats, data_name, ylog = TRUE, pvcut, rcut) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  YY <- dstats[,c("ad", "sw", "bc")]
  YY[YY < 1e-10] <- 1e-10
  ylim <- range(YY)
  nstat <- ncol(YY)
  xx <- dstats[,"rmax"]
  ylab <- expression(p[SW], p[AD], p[BC])
  xlab <- expression(R[max])
  log <- ifelse(ylog, "y", "")
  log <- rep(log, nstat)
  main <- c("Shapiro-Wilk Test",
            "Anderson-Darling Test",
            "Berkowitz Correlation Test")
  par(mfrow = c(1,3), mar = c(3.6,2.5,2,.5)+.1, oma = c(0,1,2,0))
  for(ii in 1:nstat) {
    plot(xx, YY[,ii], log = log[ii], ylim = ylim,
         xlab = "", ylab = "", pch = 16, cex = .5)
    if(!missing(pvcut)) abline(h = pvcut, lty = 2)
    if(!missing(rcut)) abline(v = rcut, lty = 2)
    title(xlab = xlab,
          #ylab = ylab[ii],
          line = 2.5)
    title(main = main[ii])
  }
  mtext(paste0("Dataset: ", data_name), side = 3,
        outer = TRUE, cex = 1, font = 2, line = .1)
  mtext("p-value", side = 2, cex = .8, outer = TRUE)
}

#--- pathwise goodness-of-fit plot ---------------------------------------------

gof_plot <- function(Xt, dT, theta, type = c("qq", "hist"), main) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  N <- nrow(Xt)-1
  tseq <- (0:N)*dT # time sequence
  type <- match.arg(type) # plot type
  Main <- main
  # residuals
  dX <- apply(Xt, 2, diff)
  Z <- fbm_resid(theta = theta, dX = dX, dT = dT)
  par(mfrow = c(2,3), mar = c(3.6, 3.6, 2, .5)+.1, oma = c(0,0,2.5,0))
  # 1D trajectories
  xlab <- expression(X[t], Y[t])
  main <- paste0("Trajectory on ", c("X", "Y"), "-axis")
  for(ii in 1:2) {
    plot(tseq[1:N], dX[,ii], main = main[ii],
         xlab = "", ylab = "", type = "l")
    title(xlab = "Time (s)", ylab = xlab[ii], line = 2.5)
  }
  # 2D trajectory
  plot(Xt[,1], Xt[,2], type = "l", xlab = "", ylab = "",
       main = "2D Trajectory")
  title(xlab = xlab[1], ylab = xlab[2], line = 2.5)
  # GoF tests
  main <- paste0("Residuals: PC", c(1, 2))
  main <- c(main, "Residuals: Combined")
  xlab <- expression(Z[1], Z[2])
  pval.ad <- c(ad_pval(Z[,1]), ad_pval(Z[,2]), ad_pval(c(Z)))
  pval.sw <- c(sw_pval(Z[,1]), sw_pval(Z[,2]), sw_pval(c(Z)))
  pval.bc <- c(bc_pval(Z[,1]), bc_pval(Z[,2]))
  pval.bc <- c(pval.bc, pbeta(min(pval.bc), 1, 2))
  for(ii in 1:3) {
    lgd <- paste0("p[", c("AD", "SW", "BC"), "]==",
                  signif(c(pval.ad[ii], pval.sw[ii], pval.bc[ii]), 2))
    if(ii <= 2) {
      zz <- Z[,ii]
    } else {
      zz <- c(Z)
    }
    if(type == "hist") {
      hist_plot(zz, xlab = xlab[ii], main = main[ii], lgd = lgd)
    } else if(type == "qq") {
      qq_plot(zz, main = main[ii], lgd = lgd)
    }
  }
  # outer margin
  mtext(Main, side = 3, outer = TRUE, cex = 1, font = 2, line = .5)
}

# helper function for histogram
hist_plot <- function(Z, xlab, main = "", lgd) {
  hist(Z, breaks = 25, freq = FALSE, main = main,
       xlab = "", ylab = "")
  curve(dnorm, add = TRUE, col = "red")
  title(xlab = xlab, ylab = "Density", line = 2.5)
  if(!missing(lgd)) legend("topright", legend = parse(text = lgd))
}

# helper function for qq plot
qq_plot <- function(Z, main = "", lgd) {
  qq <- qqnorm(Z, plot.it = FALSE)
  plot(0, type = "n", xlim = range(qq$x), ylim = range(qq$y),
       xlab = "", ylab = "", main = main)
  abline(a = 0, b = 1, lty = 2, col = "blue")
  points(qq$x, qq$y, pch = 16, cex = .5)
  title(xlab = "Theoretical Quantiles",
        ylab = "Sample Quantiles", line = 2.5)
  if(!missing(lgd)) legend("bottomright", legend = parse(text = lgd))
}

#--- msd reconstruction --------------------------------------------------------

# path reconstruction
fbm_recon <- function(theta, dX, dT, Z) {
  if(missing(Z)) Z <- fbm_resid(theta, dX, dT) # residuals
  # map residuals to corresponding normal quantiles
  N <- nrow(Z)
  q <- ncol(Z)
  nq <- if(q == 1) 1 else 3
  Zr <- matrix(qnorm(rank(Z)/(length(Z)+1)), N, ncol(Z))
  # transform to data scale
  Zr <- cholZX(Z = Zr, acf = fbm_acf(alpha = itrans_alpha(theta[1]), dT, N))
  ed <- eigen(itrans_Sigma(lambda = theta[1+q + 1:nq]))
  C <- t(ed$vec) * sqrt(ed$val)
  t(t(Zr %*% C) + theta[1+1:q] * dT)
}
