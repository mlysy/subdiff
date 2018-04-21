
# truncted estimation function for fdy
fdy_trunc <- function(dX, dT, tau) {
  # memory allocation
  N <- nrow(dX)
  qq <- ncol(dX)
  theta_hat <- rep(NA, 1)
  theta_names <- c("gamma")
  Tz <- Toeplitz(n = N)
  # profile likelihood on regular scale
  ll.prof <- function(theta) {
    alpha <- itrans_alpha(theta)
    Tz$setAcf(fdyn_acf(alpha, tau, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.prof(suff)
  }
  # calculate MLE
  theta_hat <- optimize(f = ll.prof,
                        interval = c(-5.3, -0.02), maximum = TRUE)$maximum
  names(theta_hat) <- theta_names
  # variance estimate
  V_hat <- hessian(ll.prof, x = theta_hat)
  V_hat <- subdiff:::solveV(-V_hat)
  colnames(V_hat) <- theta_names
  rownames(V_hat) <- theta_names
  list(coef = theta_hat, vcov = V_hat)
}
