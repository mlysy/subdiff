#' Base class for CSI models.
#'
#' @description Base class for CSI models.
#'
#' @importFrom stats optim optimize
#' @export
csi_model <- R6::R6Class(
  classname = "csi_model",

  private = list(
    dX_ = NULL, # internal increments.
    dt_ = NULL, # internal interobservation time.
    N_ = NULL, # internal number of increments.
    n_dims = NULL, # internal number of dimensions.
    n_drift = NULL, # internal number of drift basis coefficients.
    n_phi = NULL, # internal number of drift + acf parameters.
    omega_names = NULL, # Parameter names in the computational basis.

    Tz_ = NULL, # Toeplitz matrix object.

    #` Deep clone method.
    #`
    #` Required to create a new R6 object for `Tz_`.
    deep_clone = function(name, value) {
      switch(name,
             Tz_ = private$Tz_$clone(deep = TRUE),
             value)
    },

    #` @description Internal drift implementation.
    #`
    #` @details It seems difficult to add drift as a dynamic function with roxygen documentation.  The following workaround has the internal drift as a private member without doc, whereas the visible drift just calls the implementation.
    drift_impl = NULL,

    #` @description Calculate the profile likelihood sufficient statistics.
    #`
    get_suff = function(phi) {
      # acf and drift
      acf <- self$acf(phi, dt = private$dt_, N = private$N_)
      dr <- self$drift(phi, dt = private$dt_, N = private$N_)
      # sufficient statistics
      private$Tz_$set_acf(acf)
      lmn_suff(Y = private$dX_, X = dr, V = private$Tz_, Vtype = "acf")
    }

  ),


  active = list(

    #' @field dX Trajectory increments.  A matrix where each row is an observation and each column is a measurement dimension.  Only reallocates memory for the internal Toeplitz matrix if necessary.
    dX = function(value) {
      if(missing(value)) {
        return(private$dX_)
      } else {
        if(!is.numeric(value) || !is.matrix(value)) {
          stop("`dX` must be a numeric matrix.")
        }
        if(is.null(private$N_) || (nrow(value) != private$N_)) {
          # reallocate Toeplitz matrix
          private$N_ <- nrow(value)
          private$n_dims <- ncol(value)
          private$Tz_ <- Toeplitz$new(N = private$N_)
        }
        private$dX_ <- value
      }
    },

    #' @field dt Interobservation time (scalar).
    dt = function(value) {
      if(missing(value)) {
        return(private$dt_)
      } else {
        if(!is.numeric(value) || length(value) != 1 || value <= 0) {
          stop("`dt` must be a positive scalar.")
        }
        private$dt_ <- value
      }
    }
  ),

  public = list(
    #' @description Transform kernel parameters from regular to computational basis.
    #' @param phi Kernel parameters in the original basis.
    #' @return Kernel parameters in the computational basis.
    trans = function(phi) stop(undef_msg("trans()")),

    #' @description Transform kernel parameters from computational to regular basis.
    #'
    #' @param psi Kernel parameters in the computational basis.
    #' @return Kernel parameters in the original basis.
    itrans = function(psi) stop(undef_msg("itrans()")),

    #' @description Increment autocorrelation function.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param dt Interobservation time.
    #' @param N Number of trajectory increments.
    #' @return A vector of `N` autocorrelations.
    acf = function(phi, dt, N) stop(undef_msg("acf()")),

    #' @description Position mean square displacement function.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param t Vector of time points at which to calculate the MSD.
    #' @return A vector of MSD values the same length as `t`.
    #'
    #' @details This function can be directly supplied by the derived class, or by default is derived automatically from a high-resolution evaluation of `super$acf()`.  That is, with `dt` being up to 10x smaller than `min(diff(sort(abs(t))))`.  Thus the default method is most efficient when `t` consists of evenly spaced timepoints, and most wasteful when two time points are much closer than any others.
    #'
    #' @note The MSD of some models (e.g., `fsd` and `farma`) is not defined in continuous time, and thus depends on the interobservation time `dt`.  This is completely overlooked in `csi_model$msd()` as it is presently coded, which is potentially misleading.  Need to think carefully about what the `msd()` method should do.
    msd = function(phi, t) {
      # method 1
      # method 2
      dt <- get_dt(t)
      N <- max(abs(t)) %/% dt
      acf <- self$acf(phi, dt = dt, N = N)
      msd <- c(0, SuperGauss::acf2msd(acf))
      # pair msd with corresponding elements of tseq
      msd[abs(t) %/% dt + 1]
    },

    ## acf = function(phi, dt, N) {
    ##   if(is.null(private$acf_impl)) {
    ##     stop(undef_msg("acf()"))
    ##   } else {
    ##     private$acf_impl(phi, dt, N)
    ##   }
    ## },

    #' @description Increment autocorrelation function.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param dt Interobservation time.
    #' @param N Number of trajectory increments.
    #' @return An `N x n_drift` matrix of drift basis functions.
    drift = function(phi, dt, N) {
      if(is.null(private$drift_impl)) {
        stop(undef_msg("drift()"))
      } else {
        private$drift_impl(phi, dt, N)
      }
    },

    #' @field phi_names Vector of kernel parameter names (on the original scale).  Setting to `NULL` means there are no kernel parameters (`n_phi = 0`).
    phi_names = NA,

    #' @description Combine kernel parameters with conditional MLE of `mu` and `Sigma`.
    #'
    #' @param psi Kernel parameters in the computational basis.
    #' @return Named vector of full parameter set in the computational basis.
    get_omega = function(psi) {
      nu <- self$nu_hat(self$itrans(psi)) # nuisance terms
      setNames(c(psi, nu$mu, trans_Sigma(nu$Sigma)),
                      nm = private$omega_names)
    },

    #' @description Convert from original to computational basis.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param mu Drift coefficients.
    #' @param Sigma Scale matrix.
    #' @return Full parameter vector in the computational basis.
    trans_full = function(phi, mu, Sigma) {
      c(self$trans(phi), mu, trans_Sigma(Sigma))
    },

    #' @description Convert from computational to original basis.
    #'
    #' @param omega Vector of parameters in the computational basis.
    #' @return List with elements `phi`, `mu`, and `Sigma`.
    itrans_full = function(omega) {
      n_phi <- private$n_phi
      n_drift <- private$n_drift
      n_dims <- private$n_dims
      # number of parameters in the log-cholesky factor
      n_chol <- n_dims * (n_dims+1) / 2
      if(n_phi > 0) {
        phi <- self$itrans(omega[1:n_phi])
      } else phi <- NULL
      if(n_drift > 0) {
        mu <- matrix(omega[n_phi + 1:(n_drift*n_dims)], n_drift, n_dims)
      } else mu <- NULL
      Sigma <- itrans_Sigma(omega[n_phi + n_drift*n_dims + 1:n_chol])
      list(phi = phi, mu = mu, Sigma = Sigma)
    },

    #' @description Evaluate the negative profile loglikelihood.
    #'
    #' @param psi Kernel parameters in the computational basis.
    #' @return Value of the negative profile loglikelihood (scalar).
    nlp = function(psi) {
      phi <- self$itrans(psi) # convert psi to original scale
      suff <- private$get_suff(phi) # sufficient statistics
      -lmn_prof(suff)
    },

    #' @description Conditional MLE of nuisance parameters.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @return A list with elements `mu` and `Sigma`.
    nu_hat = function(phi) {
      suff <- private$get_suff(phi) # sufficient statistics
      list(mu = suff$Bhat, Sigma = suff$S/suff$n)
    },

    #' @description Evaluate the loglikelihood function.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param mu Drift coefficients.
    #' @param Sigma Scale matrix.
    #' @return Value of the loglikelihood (scalar).
    loglik = function(phi, mu, Sigma) {
      suff <- private$get_suff(phi) # sufficient statistics
      lmn_loglik(Beta = mu, Sigma = Sigma, suff = suff)
    },

    #' @description Calculate the observed Fisher information matrix.
    #'
    #' @param omega Vector of length `n_omega` of parameters in the computational basis.
    #' @return The `n_omega x n_omega` observed Fisher information matrix.
    fisher = function(omega) {
      numDeriv::hessian(x = omega, func = function(omega) {
        # convert omega to phi, mu, Sigma
        theta <- self$itrans_full(omega)
        -self$loglik(theta$phi, theta$mu, theta$Sigma)
      })
    },

    #' @description Convert a Fisher information matrix to a variance matrix.
    #'
    #' @param fi Fisher information matrix of size `n_omega x n_omega` with parameters in the computational basis.
    #' @return Variance matrix of size `n_omega x n_omega` with parameters in the computational basis.
    get_vcov = function(fi) {
      # invert to get variance estimate
      var_hat <- chol2inv(chol(fi))
      colnames(var_hat) <- rownames(var_hat) <- private$omega_names
      var_hat
    },

    #' @description Calculate the maximum likelihood parameter values.
    #'
    #' @param psi0 Vector of kernel parameter values (on the computational scale) to initialize the optimization if `n_phi > 1`.  If `n_phi == 1`, a vector of length 2 giving the range in which to perform the optimum search.
    #' @param vcov Whether to also calculate the MLE variance estimate.
    #' @param ... Additional arguments to [stats::optim()] or [stats::optimize()] for `n_phi == 1`.
    #'
    #' @return A list with elements `coef` and optionally `vcov` containing the MLE in the computational basis and its variance estimate.
    fit = function(psi0, vcov = TRUE, ...) {
      # calculate MLE
      if(private$n_phi > 1) {
        opt <- optim(par = psi0,
                     fn = self$nlp, ...)
        if(opt$convergence != 0) warning("`optim()` did not converge.")
        psi_hat <- opt$par
      } else {
        psi_hat <- optimize(f = self$nlp, interval = psi0, ...)$minimum
      }
      omega_hat <- self$get_omega(psi_hat)
      ## phi_hat <- self$itrans(psi_hat)
      ## # nuisance terms
      ## nu <- self$nu_hat(phi_hat)
      ## omega_hat <- c(psi_hat, nu$mu, trans_Sigma(nu$Sigma))
      ## names(omega_hat) <- private$omega_names
      if(vcov) {
        # fisher information
        fi <- self$fisher(omega_hat)
        # invert to get variance estimate
        var_hat <- self$get_vcov(fi)
        ## var_hat <- chol2inv(chol(fi))
        ## colnames(var_hat) <- rownames(var_hat) <- private$omega_names
        out <- list(coef = omega_hat, vcov = var_hat)
      } else out <- omega_hat
      out
    },

    #' @description Calculate the model residuals.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param mu Drift coefficients.
    #' @param Sigma Scale matrix.
    #'
    #' @return A matrix of residuals the same size as `dX` as calculated with [csi_resid()], upon using the model's `drift()` and `acf()` specifications.
    resid = function(phi, mu, Sigma) {
      dr <- self$drift(phi, dt = self$dt, N = private$N_) %*% mu
      ac <- self$acf(phi, dt = self$dt, N = private$N_)
      csi_resid(self$dX, drift = dr, acf = ac, Sigma = Sigma)
    },

    #' @description Model object constructor.
    #'
    #' @param dX Trajectory increments.
    #' @param dt Interobservation time.
    #' @param drift Drift specification.  Either one of the strings "none", "linear", "quadratric", or a function with signature `function(phi, dt, N)`.
    #' @param n_drift Integer number of drift terms.  Ignored if `drift` is one of the default strings.  Required otherwise.
    #'
    #' @details
    #' - Is it worth checking whether model object is valid at construction time?
    #' - `n_phi` is automatically determined from `phi_names`.  But this means the latter must be set by `derived$initialize()` before `super$initialize()` is called.  Otherwise an error is thrown.
    initialize = function(dX, dt, drift = "linear", n_drift) {
      if(!missing(dX)) self$dX <- dX
      self$dt <- dt
      # set drift
      check_drift(drift) # validate
      if(drift == "linear") {
        private$drift_impl <- drift_linear
        private$n_drift <- 1
      } else if(drift == "none") {
        private$drift_impl <- drift_none
        private$n_drift <- 0
      } else if(drift == "quadratic") {
        private$drift_impl <- drift_quadratic
        private$n_drift <- 2
      } else {
        # custom drift
        private$drift_impl <- drift
        private$n_drift <- n_drift
      }
      # kernel parameters
      if(anyNA(self$phi_names)) stop(undef_msg("phi_names")) # must have already been set
      private$n_phi <- length(self$phi_names)
      # computational parameter names
      private$omega_names <- get_omega_names(n_dims = private$n_dims,
                                             n_drift = private$n_drift,
                                             n_phi = private$n_phi)
    }
  )
)

#--- helper functions ----------------------------------------------------------

#' Check that drift arguments are correctly specified.
#'
#' @param drift Drift function.
#' @noRd
check_drift <- function(drift) {
  if(length(drift) != 1 || (!is.character(drift) &&
                            !is.function(drift))) {
    stop("Incorrect drift specification.  Please see `csi_model` for details.")
  }
  if(is.character(drift) && !(drift %in% c("linear", "none", "quadratic"))) {
    stop('Prespecified drift must be one of "linear", "none", or "quadratic".')
  }
  if(is.function(drift) && !identical(methods::formalArgs(drift), c("phi", "dt", "N"))) {
    stop("`drift()` must have argument signature `phi`, `dt`, `N`.")
  }
}

#' Generic error message if a concrete specification of an abstract member/method has not been set.
#'
#' @param name Name of abstract member/method.
#' @noRd
undef_msg <- function(name) {
  paste0("`self$", name, "` has not been defined yet.")
}

#' Set computational basis parameter names.
#'
#' @param n_dims Number of trajectory dimensions.
#' @param n_drift Number of drift dimensions.
#' @param n_phi Number of act parameters.
#'
#' @details The computational basis parameters are named in the following order:
#' - `psi`: Numbered from 1 to `n_phi`.
#' - `mu`: Double indexing from the matrix of size `n_drift x n_dims` in column-major order.
#' - `lambda`: Double indexing from upper triangular matrix of size `n_dims x n_dims`.
#' @noRd
get_omega_names <- function(n_dims, n_drift, n_phi) {
  if(n_phi > 0) {
    psi_names <- paste0("psi", 1:n_phi)
  } else psi_names <- NULL
  if(n_drift > 0) {
    mu_names <- paste0("mu", apply(expand.grid(1:n_drift, 1:n_dims), 1,
                                   paste0, collapse = ""))
  } else mu_names <- NULL
  lambda_names <- paste0("lambda", apply(expand.grid(1:n_dims, 1:n_dims), 1,
                                         paste0, collapse = ""))
  lambda_names <- matrix(lambda_names, n_dims, n_dims)
  lambda_names <- lambda_names[upper.tri(lambda_names, diag = TRUE)]
  c(psi_names, mu_names, lambda_names)
}

#' Pick a sampling frequency which is sufficiently close to each element of `t`.
#'
#' @param t Vector of positive time points.
#' @param tol Tolerance (see 'Details').
#' @param n Number of subdivisions between `tol` and `min(t)`
#' @details `dt` will be such that `max(t modulo dt) < tol`.
#' @noRd
get_dt <- function(t, tol, n = 10) {
  tseq <- sort(abs(t))
  tseq <- tseq[tseq > .Machine$double.eps]
  min_diff <- min(diff(c(0, tseq)))
  if(missing(tol)) {
    # pick a minimum tolerance
    # fixme: what if tseq[1] == 0?
    tol <- min_diff/10
  }
  dt_seq <- sort(seq(min_diff, tol, len = n))
  dt_err <- sapply(dt_seq, function(dt) max(tseq %% dt))
  dt_seq[which.max(dt_err < tol)]
}
