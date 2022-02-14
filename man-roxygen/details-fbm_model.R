#' @details Fractional Brownian motion `X_t` is the only (zero-mean) continuous stationary increments (CSI) Gaussian process having mean square displacement (MSD) function given by the power law
#' ```
#' E[(X_t-X_0)^2] = t^alpha,
#' ```
#' with subdiffusion exponent `0 < alpha < 2`.  The resulting MSD as a function of `alpha` is what gets passed to the [csi_model] base class to construct the `fbm_model` derived class.
