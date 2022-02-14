#' Base class for CSI models.
#'
#' @description Base class for CSI models.
#'
#' @details Let `Xt` denote an `N x d` matrix of particle positions recorded at intervals of `dt`, where `d = 1,2,3` is the number of recorded dimensions.  A CSI model for `Xt` is of the form
#' ```
#' Xt = R(phi) mu + Sigma^{1/2} Z,
#' ```
#' where:
#'
#' - `R(phi)` is an `N x p` matrix of drift terms, possibly dependent on a parameter `phi`.  In most cases though, `R(phi) = t((1:N-1) * dt)` is used to model linear drift.
#' - `mu` is a `p x d` matrix of drift parameters.  For linear drift, it represents the drift velocity in each of the `d` dimensions.
#' - `Sigma` is a `d x d` variance matrix.
#' - `Z` is an `N x d` matrix where each column is an iid realization of `N` values of a Gaussian continuous stationary increments (CSI) process.  Such processes are completely determined by the mean square displacement (MSD) function
#'     ```
#'     MSD_Z(h) = E[ (Z_{i+h,j}, Z_{i,j})^2 ],
#'     ```
#'     where the MSD function also depends on the parameter `phi`.
#'
#' The `csi_model` class is a base class for CSI models which requires the user to specify the drift and msd functions, based on which it provides generic methods for parameter inference, simulation, etc.  For further details about the CSI model see `vignette("subdiff")`.
#'
#' @importFrom stats optim optimize
#' @export
csi_model <- R6::R6Class(
  classname = "csi_model",

  private = list(
    Xt_ = NULL, # internal positions.
    dX_ = NULL, # internal increments.
    dt_ = NULL, # internal interobservation time.
    N_ = NULL, # internal number of increments.
    n_dims = NULL, # internal number of dimensions.
    n_drift = NULL, # internal number of drift basis coefficients.
    n_phi = NULL, # internal number of drift + acf parameters.
    omega_names = NULL, # Full parameter names in the computational basis.
    psi_names = NULL, # Kernel parameter names in the computational basis.

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

    #' @field Xt Particle trajectory.  A matrix where row `n` is the position of the particle at time `t = n * dt` and each column is a measurement dimension.  Only reallocates memory for the internal Toeplitz matrix if necessary.
    Xt = function(value) {
      if(missing(value)) {
        return(private$Xt_)
      } else {
        value <- check_Xt(value)
        private$Xt_ <- value
        private$n_dims <- ncol(value)
        if(is.null(private$N_) || (nrow(value) != private$N_+1)) {
          # reallocate Toeplitz matrix
          private$N_ <- nrow(value) - 1
          private$Tz_ <- Toeplitz$new(N = private$N_)
        }
        # set increments
        private$dX_ <- matrix(apply(value, 2, diff), ncol = private$n_dims)
      }
    },

    ## #' @field dX Trajectory increments.  A matrix where each row is an observation and each column is a measurement dimension.  Only reallocates memory for the internal Toeplitz matrix if necessary.
    ## dX = function(value) {
    ##   if(missing(value)) {
    ##     return(private$dX_)
    ##   } else {
    ##     if(!is.numeric(value) || !is.matrix(value)) {
    ##       stop("`dX` must be a numeric matrix.")
    ##     }
    ##     if(is.null(private$N_) || (nrow(value) != private$N_)) {
    ##       # reallocate Toeplitz matrix
    ##       private$N_ <- nrow(value)
    ##       private$n_dims <- ncol(value)
    ##       private$Tz_ <- Toeplitz$new(N = private$N_)
    ##     }
    ##     private$dX_ <- value
    ##   }
    ## },

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
    #' @return A vector of *unscaled* MSD values the same length as `t`.  That is, for the drift-subtracted process `X_sub(t) = X(t) - drift(t) * mu`, the method returns `eta(t)` where
    #'
    #' ```
    #' MSD_{X_sub}(t) = trace(Sigma) * eta(t).
    #'
    #' @details This method can be directly supplied by the derived class.  Otherwise, it uses [SuperGauss::acf2msd()] to calculate the MSD from `self$acf()` at intervals of `self$dt`, and interpolates linearly between these timepoints at the desired values in `t`.
    #'
    #' The `self$msd()` method is defined this way because the MSD of some CSI models (e.g., `fsd` and `farma`) is not defined in continuous time, and thus intrinsically depends on the interobservation time `dt`.  Thus, the default `self$msd()` method throws an error if `self$dt` has not yet been set.
    msd = function(phi, t) {
      # method 1
      N <- ceiling(max(abs(t)) / self$dt)
      acf <- self$acf(phi, dt = self$dt, N = N)
      tseq <- self$dt * 1:N
      msd <- SuperGauss::acf2msd(acf)
      approx(x = c(-rev(tseq), 0, tseq),
             y = c(rev(msd), 0, msd),
             xout = t)$y
      ## # method 2
      ## dt <- get_dt(t)
      ## N <- max(abs(t)) %/% dt
      ## acf <- self$acf(phi, dt = dt, N = N)
      ## msd <- c(0, SuperGauss::acf2msd(acf))
      ## # pair msd with corresponding elements of tseq
      ## msd[abs(t) %/% dt + 1]
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
      setNames(c(self$trans(phi), mu, trans_Sigma(Sigma)),
               nm = private$omega_names)
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
      fi <- numDeriv::hessian(x = omega, func = function(omega) {
        # convert omega to phi, mu, Sigma
        theta <- self$itrans_full(omega)
        -self$loglik(theta$phi, theta$mu, theta$Sigma)
      })
      colnames(fi) <- rownames(fi) <- private$omega_names
      fi
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
      csi_resid(private$dX_, drift = dr, acf = ac, Sigma = Sigma)
    },

    #' @description Simulate trajectories from the model.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param mu Drift coefficients.
    #' @param Sigma Scale matrix.
    #' @param nsim Number of trajectories to simulate.
    #' @param fft,nkeep,tol Optional arguments to [SuperGauss::rnormtz()].
    #'
    #' @return A matrix of size `dim(Xt)` or an array of size `nrow(Xt) x ncol(dX) x nsim` array when `nsim > 1` of simulated trajectories, as calculated with [csi_sim()], using the model's `drift()` and `acf()` specifications.
    sim = function(phi, mu, Sigma, nsim = 1, fft = TRUE, nkeep, tol = 1e-6) {
      dr <- self$drift(phi, dt = self$dt, N = private$N_) %*% mu
      ac <- self$acf(phi, dt = self$dt, N = private$N_)
      csi_sim(drift = dr, acf = ac, Sigma = Sigma, X0 = self$Xt[1,],
              nsim = nsim, fft = fft, nkeep = nkeep, tol = tol)
    },

    #' @description Model object constructor.
    #'
    #' @param Xt Matrix of particle positions.
    #' @param dt Interobservation time.
    #' @param drift Drift specification.  Either one of the strings "none", "linear", "quadratric", or a function with signature `function(phi, dt, N)`.
    #' @param n_drift Integer number of drift terms.  Ignored if `drift` is one of the default strings.  Required otherwise.
    #'
    #' @details
    #' The value of `n_phi` is automatically determined from `phi_names`.  But this means the latter must be set by `derived$initialize()` before `super$initialize()` is called.  Otherwise an error is thrown.
    #'
    #' The constructor can be called without `Xt` for accessing methods which don't require it, e.g.,
    #' ```
    #' derived$new(dt = dt)$acf(phi, dt, N)
    #' ```
    #'
    #' *Development Notes*
    #'
    #' - Should a validator be (optionally) called after the constructor?
    #' - Should `psi` necessarily contain the origin, to facilitate the validator?
    initialize = function(Xt, dt, drift = "linear", n_drift) {
      if(!missing(Xt)) self$Xt <- Xt
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
      private$psi_names <- private$omega_names[1:private$n_phi]
    }
  )
)

#--- helper functions ----------------------------------------------------------

#' Check that `Xt` is a numeric matrix.
#'
#' If it's a data frame, convert to numeric matrix.
#'
#' @template args-Xt
#' @noRd
check_Xt <- function(Xt) {
  if(is.list(Xt)) Xt <- as.matrix(Xt)
  if(!is.numeric(Xt) || !is.matrix(Xt)) {
    stop("`Xt` must be a numeric matrix.")
  }
  Xt
}

#' Check that drift arguments are correctly specified.
#'
#' @param drift Drift function.
#' @noRd
check_drift <- function(drift) {
  if(length(drift) != 1 || (!is.character(drift) &&
                            !is.function(drift))) {
    stop("Incorrect drift specification.  Please see [csi_model$drift()] for details.")
  }
  if(is.character(drift) && !(drift %in% c("linear", "none", "quadratic"))) {
    stop('Prespecified drift must be one of "linear", "none", or "quadratic".')
  }
  if(is.function(drift) && !identical(methods::formalArgs(drift), c("phi", "dt", "N"))) {
    stop("`drift()` must have argument signature `function(phi, dt, N)`.")
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
