#' Class definition for the fSD model.
#'
#' @name fsd_model
#'
#' @details Constructor function that generates a list of fsd model for purpose of fitting.
#'
#' @example examples/fsd_model.R
#'
NULL

#' @rdname fsd_model
#' @export
fsd_model <- R6::R6Class(
  classname = "fsd_model",
  inherit = csi_model,

  public = list(
    #' @field phi_names Kernel parameter names.  A subset of `(alpha, tau, sigma2)` (see [csi_model] and `fsd_model$initialize()`).
    phi_names = c("alpha", "tau", "sigma2"),

    #' @description Increment autocorrelation function.
    #'
    #' @param phi,dt,N See [csi_model].
    acf = function(phi, dt, N) {
      fsd_acf(alpha = phi[1], tau = phi[2], sigma2 = phi[3], dt = dt, N = N)
    },

    #' @description Transform kernel parameters from regular to computational basis.
    #'
    #' @param phi See [csi_model].
    #' @details The transformation function is
    #' ```
    #' psi = (psi1 = logit(alpha/2), psi2 = logit(tau), psi3 = log(sigma2)).
    #' ```
    trans = function(phi) {
      setNames(c(logit(phi[1], min = 0, max = 2), # alpha
                 logit(phi[2], min = 0, max = 1), # tau
                 log(phi[3])), # sigma2
               nm = paste0("psi", 1:private$n_phi))
    },

    #' @description Transform kernel parameters from computational to regular basis.
    #'
    #' @param psi See [csi_model].
    itrans = function(psi) {
      setNames(c(ilogit(psi[1], min = 0, max = 2), # alpha
                 ilogit(psi[2], min = 0, max = 1), # tau
                 exp(psi[3])), # sigma2
               nm = self$phi_names)
    },

    #' @description Transform parameters from computational basis to subdiffusion parameters.
    #'
    #' @param omega See [csi_model].
    #' @return Vector with named elements `alpha` and `logD`.
    get_subdiff = function(omega) {
      theta <- self$itrans_full(omega) # convert to inferential basis
      # extract alpha and logD
      setNames(c(theta$phi["alpha"],
                 log(mean(diag(theta$Sigma)))),
               nm = c("alpha", "logD"))
    }
  )

)

## fsd_model <- function() {
##   # parameter name
##   fsd_theta <- c("alpha", "tau", "sigma2")

##   # ACF function
##   .fsd_acf <- function(theta, dt, N) {
##     fsd_acf(theta[1], theta[2], theta[3], dt, N)
##   }

##   # transformation, from original scale (theta) to unrestricted scale (gamma)
##   fsd_trans <- function(theta) {
##     gamma <- theta
##     # alpha
##     gamma[1] <- logit(theta[1], min = 0, max = 2)
##     # tau
##     gamma[2] <- logit(theta[2], min = 0, max = 1)
##     # sigma2
##     gamma[3] <- log(theta[3])
##     gamma
##   }

##   # inverse transformation, from unrestricted scale (gamma) to original scale (theta)
##   fsd_itrans <- function(gamma) {
##     theta <- gamma
##     # alpha
##     theta[1] <- ilogit(gamma[1], min = 0, max = 2)
##     # tau
##     theta[2] <- ilogit(gamma[2], min = 0, max = 1)
##     # sigma2
##     theta[3] <- exp(gamma[3])
##     theta
##   }

##   # penalty function on transformed scale
##   fsd_penalty <- function(gamma) {
##     log1pe(gamma[2]) + log1pe(-gamma[2]) - gamma[3]
##   }

##   # plug-in values
##   # fsd_plugin <- names_plugin <- NULL
##   # if(!missing(tau)) {
##   #   fsd_plugin <- c(fsd_plugin, tau)
##   #   names_plugin <- c(names_plugin, "tau")
##   # }
##   # if(!missing(sigma2)) {
##   #   fsd_plugin <- c(fsd_plugin, sigma2)
##   #   names_plugin <- c(names_plugin, "sigma2")
##   # }
##   # names(fsd_plugin) <- names_plugin

##   model <- list(theta_names = fsd_theta,
##                 acf = .fsd_acf,
##                 theta_trans = fsd_trans,
##                 theta_itrans = fsd_itrans,
##                 penalty = fsd_penalty)
##   class(model) <- "csi_model"
##   model
## }
