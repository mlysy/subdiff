
source(system.file("proj", "plot-functions.R", package = "subdiff"))

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

#-------------------------------------------------------------------------------

# simulate data from fbm model using cholesky method
fbm_sim <- function(theta, dT, N, Z) {
  alpha <- itrans_alpha(theta[1])
  mu <- theta[2:3]
  Sigma <- itrans_Sigma(theta[4:6])
  acf <- fbm_acf(alpha, dT, N-1)
  if(missing(Z)) {
    dX <- rSnorm(2, acf) %*% chol(Sigma)
  } else {
    if(any(dim(Z) != c(N-1, 2))) stop("Z has wrong dimensions.")
    dX <- cholZX(Z = Z, acf = acf)
    ed <- eigen(Sigma)
    C <- t(ed$vec) * sqrt(ed$val)
    dX <- dX %*% C
    ## dX <- t(t(dX %*% C) + mu * dT)
  }
  apply(rbind(0, t(t(dX) + mu * dT)), 2, cumsum)
}

#' fBM path reconstruction.
#'
#' Deterministic mapping of a given path to an idealized version under the fBM+Drift model.
#'
#' @param theta Parameters of fBM + drift model.
#' @param dX Vector or 2-column matrix of increments.
#' @param dT Interobservation time (scalar).
#' @param Z Optional white noise residuals having the same dimensions as \code{dX}.  Defaults to model residuals calculated with \code{theta} (see Details).
#' @return Reconstructed trajectory, having the same dimensions as \code{dX} (or \code{Z}).
#' @details The incremements of the fBM + drift model can be written as
#' \preformatted{
#' dX = V_alpha^{1/2} Z Sigma^{-1/2} + mu dT,    Z_{ij} ~iid N(0,1),
#' }
#' where \code{theta = (alpha, mu, Sigma)} are the parameters of the model, and we can take any square root of the covariance matrices \code{Sigma} and \code{V_alpha} corresponding to spatial and temporal dependence.  The choice of eigendecomposition for the former, and Cholesky for the latter is made here.
#'
#' Given a set of residuals \code{Z}, the idealized residuals correspond to the expected values of the \code{N(0,1)} order statistics corresponding to the original quantiles' ranks.
fbm_recon <- function(theta, dX, dT, Z) {
  if(missing(Z)) Z <- fbm_resid(theta, dX, dT) # residuals
  # map residuals to corresponding normal quantiles
  N <- nrow(Z)
  qq <- ncol(Z)
  nq <- if(qq == 1) 1 else 3
  Zr <- matrix(qnorm(rank(Z)/(length(Z)+1)), N, ncol(Z))
  # transform to data scale
  Zr <- cholZX(Z = Zr, acf = fbm_acf(alpha = itrans_alpha(theta[1]), dT, N))
  ed <- eigen(itrans_Sigma(lambda = theta[1+qq + 1:nq]))
  C <- t(ed$vec) * sqrt(ed$val)
  t(t(Zr %*% C) + theta[1+1:qq] * dT)
}

