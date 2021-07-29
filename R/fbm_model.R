#' Class definition for the fBM model.
#'
#' @name fbm_model
#'
#' @details Fractional Brownian motion `X_t` is the only (zero-mean) continuous stationary increments (CSI) Gaussian process having mean square displacement (MSD) function given by the power law
#' ```
#' E[(X_t-X_0)^2] = t^alpha,
#' ```
#' with subdiffusion exponent `0 < alpha < 2`.
#'
#' @example examples/fbm_sim.R
#' @example examples/fbm_model.R
NULL

#' @rdname fbm_model
#' @export
fbm_model <- R6::R6Class(
  classname = "fbm_model",
  inherit = csi_model,

  public = list(
    #' @field phi_names Kernel parameter names. In this case, the character string `alpha`.  See [csi_model].
    phi_names = "alpha",
    #' @description Increment autocorrelation function.
    #'
    #' @param phi,dt,N See [csi_model].
    acf = function(phi, dt, N) fbm_acf(alpha = phi, dt = dt, N = N),
    #' @description Transform kernel parameters from regular to computational basis.
    #'
    #' @param phi See [csi_model].
    trans = function(phi) setNames(logit(phi, 0, 2), nm = "psi1"),
    #' @description Transform kernel parameters from computational to regular basis.
    #'
    #' @param psi See [csi_model].
    itrans = function(psi) setNames(ilogit(psi, 0, 2), nm = self$phi_names),

    #' @description Transform parameters from computational basis to subdiffusion parameters.
    #'
    #' @param omega See [csi_model].
    #' @return Vector with named elements `alpha` and `logD`.
    get_subdiff = function(omega) {
      theta <- self$itrans_full(omega) # convert to inferential basis
      # extract alpha and logD
      setNames(c(theta$phi["alpha"],
                 log(mean(diag(theta$Sigma))/2)),
               nm = c("alpha", "logD"))
    }
  )
)

## fbm_model$set("public", "phi_names", "alpha", overwrite = TRUE)
## fbm_model$set("public", "acf",
##               function(phi, dt, N) fbm_acf(alpha = phi, dt = dt, N = N),
##               overwrite = TRUE)
## fbm_model$set("public", "trans",
##               function(phi) logit(phi, 0, 2),
##               overwrite = TRUE)
## fbm_model$set("public", "itrans",
##               function(psi) ilogit(psi, 0, 2),
##               overwrite = TRUE)

## fbm_model <- function() {
##   fbm_theta <- "alpha"
##   fbm_trans <- function(alpha) logit(alpha, min = 0, max = 2)
##   fbm_itrans <- function(theta) ilogit(theta, min = 0, max = 2)
##   fbm_penalty <- function(gamma) 0
##   fbm_penalty <- function(gamma) 0

##   model <- list(theta_names = fbm_theta,
##                 acf = fbm_acf,
##                 theta_trans = fbm_trans,
##                 theta_itrans = fbm_itrans,
##                 penalty = fbm_penalty)
##   class(model) <- "csi_model"
##   model
## }
