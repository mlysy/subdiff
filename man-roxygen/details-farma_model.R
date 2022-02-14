#' @details Let `X_n` denote the true position of an fBM process at time `t = n dt`.  The farma(p,q) error model describes the measured position `Y_n` at time `t = n dt` as
#' ```
#' Y_n = sum_{i=1}^p phi_i Y_{n-i} + sum_{j=0}^q rho_j X_{n-j},
#' ```
#' with `rho_0 = 1 - sum_{i=1}^p phi_i - sum_{j=1}^q rho_j`.  The resulting MSD as a function of `phi` and `rho` is what gets passed to the [csi_model] base class to construct the `farma_model` derived class.
