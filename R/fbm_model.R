#' Constructor of fBM model object.
#'
#' @name fbm_model
#' @details Constructor function that generates a list of fBM model for purpose of fitting.
#' @template ret-csi_class
#'
#' @example examples/fbm_model.R
NULL

#' @rdname fbm_model
#' @export
fbm_model <- R6::R6Class(
  classname = "fbm_model",
  inherit = csi_model,

  public = list(
    #' @field phi_names See [csi_model].
    phi_names = "alpha",
    #' @description Increment autocorrelation function.
    #'
    #' @param phi,dt,N See [csi_model].
    acf = function(phi, dt, N) fbm_acf(alpha = phi, dt = dt, N = N),
    #' @description Transform kernel parameters from regular to computational basis.
    #'
    #' @param phi See [csi_model].
    trans = function(phi) logit(phi, 0, 2),
    #' @description Transform kernel parameters from computational to regular basis.
    #'
    #' @param psi See [csi_model].
    itrans = function(psi) ilogit(psi, 0, 2)
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
