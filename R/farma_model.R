#' Constructor for farma model object.
#'
#' @name farma_model
#'
#' @details Constructor function that generates a list of farma(p,q) model for purpose of fitting.
#' @template ret-csi_class
#'
#' @example examples/farma_sim.R
#' @example examples/farma_model.R
NULL

#' @rdname farma_model
#' @export
farma_model <- R6::R6Class(
  classname = "farma_model",
  inherit = csi_model,

  private = list(
    p_ = NULL, # number of AR terms.
    q_ = NULL, # number of MA terms.
    m_ = NULL, # number of terms in MA approximation.

    #` Obtain the `p` and `q` parameters.
    #`
    #` Gives an error if these are not formatted correctly.
    get_pq = function(order) {
      if(length(order) != 2) stop("order must be a vector of length 2.")
      p <- order[1]
      q <- order[2]
      ## if(missing(p)) p <- 0
      if((p - as.integer(p) != 0) && p < 0) {
        stop("AR order term `p = order[1]` must be a nonnegative integer.")
      }
      ## if(missing(q)) q <- 0
      if((q - as.integer(q) != 0) && q < 0) {
        stop("MA order term `q = order[2]` must be a nonnegative integer.")
      }
      setNames(c(p, q), nm = c("p", "q"))
    },

    #` @description Internal acf implementation.
    #`
    #` @details Workaround for dynamic members with roxygen documentation.  See [`csi_model$drift()`][csi_model].
    acf_impl = NULL
  ),

  public = list(
    #' @field phi_names Kernel parameter names.  The character vector `(alpha, phi_1, ..., phi_p, rho_1, ..., rho_q)`.
    phi_names = NA,

    #' @description Increment autocorrelation function.
    #'
    #' @param phi,dt,N See [csi_model].
    acf = function(phi, dt, N) {
      if(is.null(private$acf_impl)) {
        stop(undef_msg("acf()"))
      } else {
        private$acf_impl(phi, dt, N)
      }
    },

    #' @description Transform kernel parameters from regular to computational basis.
    #'
    #' @param phi See [csi_model].
    trans = function(phi) {
      psi <- phi
      psi[1] <- logit(phi[1], min = 0, max = 2)
      psi[-1] <- logit(phi[-1], min = -1, max = 1)
      setNames(psi, nm = paste0("psi", 1:private$n_phi))
    },

    #' @description Transform kernel parameters from computational to regular basis.
    #'
    #' @param psi See [csi_model].
    itrans = function(psi) {
      phi <- psi
      phi[1] <- ilogit(psi[1], min = 0, max = 2)
      phi[-1] <- ilogit(psi[-1], min = -1, max = 1)
      setNames(phi, nm = self$phi_names)
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
    },

    #' @description Class constructor.
    #'
    #' @param Xt,dt,drift,n_drift See [csi_model].
    #' @param order Vector of two nonnegative integers specifying the number of autoregressive and moving-average terms, respectively.
    #' @param m Order of the moving-average approximation (see [farma_acf()]).
    initialize = function(Xt, dt, order = c(0, 0), m = 50, drift = "linear", n_drift) {
      # determine p and q
      pq <- private$get_pq(order)
      p <- pq["p"]
      q <- pq["q"]
      private$p_ <- p
      private$q_ <- q
      private$m_ <- m
      # set parameter names
      phi_names <- "alpha"
      if(p > 0) phi_names <- c(phi_names, paste0("phi", 1:p))
      if(q > 0) phi_names <- c(phi_names, paste0("rho", 1:q))
      self$phi_names <- phi_names
      # create the acf function
      private$acf_impl <- function(phi, dt, N) {
        alpha <- phi[1]
        ar_coef <- if(private$p_ == 0) numeric() else phi[1+1:p]
        ma_coef <- if(private$q_ == 0) numeric() else phi[1+p+1:q]
        farma_acf(alpha = alpha, phi = ar_coef, rho = ma_coef,
                  dt = dt, N = N, m = private$m_)
      }
      super$initialize(Xt = Xt, dt = dt, drift = drift, n_drift = n_drift)
      ## FIXME: modify the drift term
    }
  )
)


## farma_model <- function(p, q) {
##   if(!p & !q) {
##     # depreciated to fBM model
##     return(fbm_model())
##   }

##   farma_theta <- "alpha"
##   if(p) farma_theta <- c(farma_theta, paste0("phi", 1:p))
##   if(q) farma_theta <- c(farma_theta, paste0("rho", 1:q))

##   .farma_acf <- function(theta, dt, N, m = 30) {
##     if(!p) {
##       farma_acf(theta[1], numeric(), theta[1+1:q], dt, N, m)
##     } else {
##       farma_acf(theta[1], theta[1+1:p], theta[1+p+1:q], dt, N, m)
##     }
##   }

##   farma_trans <- function(theta) {
##     gamma <- theta
##     # alpha
##     gamma[1] <- logit(theta[1], min = 0, max = 2)
##     # phi and rho
##     gamma[-1] <- logit(theta[-1], min = -1, max = 1)
##     gamma
##   }

##   farma_itrans <- function(gamma) {
##     theta <- gamma
##     # alpha
##     theta[1] <- ilogit(gamma[1], min = 0, max = 2)
##     # phi and rho
##     theta[-1] <- ilogit(gamma[-1], min = -1, max = 1)
##     theta
##   }

##   farma_penalty <- function(gamma) 0

##   model <- list(theta_names = farma_theta,
##                 acf = .farma_acf,
##                 theta_trans = farma_trans,
##                 theta_itrans = farma_itrans,
##                 penalty = farma_penalty)
##   class(model) <- "csi_model"
##   model
## }
