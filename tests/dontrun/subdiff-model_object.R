#' ---
#' title: "Model Objects in **subdiff**"
#' author: "Martin Lysy"
#' date: "`r Sys.Date()`"
#' output:
#'   rmarkdown::html_document:
#'     toc: yes
#' ---
#'
#' \newcommand{\bm}[1]{\boldsymbol{#1}}
#' \newcommand{\RR}{\bm{R}}
#' \newcommand{\XX}{\bm{X}}
#' \newcommand{\YY}{\bm{Y}}
#' \newcommand{\ZZ}{\bm{Z}}
#' \newcommand{\dt}{\Delta t}
#' \newcommand{\N}{\mathcal{N}}
#' \newcommand{\iid}{\stackrel{\mathrm{iid}}{\sim}}
#' \newcommand{\msd}{\operatorname{MSD}}
#' \newcommand{\var}{\operatorname{var}}
#' \newcommand{\mmu}{\bm{\mu}}
#' \newcommand{\eet}{\bm{\eta}}
#' \newcommand{\tth}{\bm{\theta}}
#' \newcommand{\lla}{\bm{\lambda}}
#' \newcommand{\pps}{\bm{\psi}}
#' \newcommand{\pph}{\bm{\phi}}
#' \newcommand{\SSi}{\bm{\Sigma}}

#+ setup, include = FALSE
knitr::opts_chunk$set(eval = FALSE)

#' # Naming Conventions
#'
#' - Model: $\XX(t) = \RR(t \mid \pph) \mmu + \SSi^{1/2} \ZZ(t), \qquad \msd_Z(t) = \eta(t \mid \pph)$.
#'
#' - Original basis: $\tth = (\pph, \mmu, \SSi)$.
#'
#' - Computational basis: $\pps = \bm{g}(\pph)$, $\eet = (\pps, \mmu, \lla)$.
#'
#' # Usage
#'
#' Model objects will inherit from a base class `csi_model` which provides a number of default methods outlined below.
#'
#'
#' ## Constructor
#'
#' - Performs memory allocation.
#' - Can optionally override the drift which defaults to "linear".
#' - Should the data be provided on the increment or regular scale?
#'
#' Examples:

#+ constructor
# direct inheritance of base class constructor. drift can be skipped.
mod <- fbm_model(dX, dt, drift = "linear")

# custom derived class constructor with additional arguments
mod <- farma_model(dX, dt, p = 3, q = 2)
mod <- fsd_model(dX, dt, tau = 0) # savin-doyle model, pure static error

#' ## Internal Variables
#'
#' A number of members and methods in `csi_model` will need to be set by the derived model class:
#'
#' - `n_drift`: The number of basis functions in the drift term.  E.g., for linear drift this is `n_drift = 1` and for quadratic drift this is `n_drift = 2`.
#' - `phi_names`: The names of the MSD and drift parameters in the original basis.
#' - `n_phi`: The dimension of $\pph$.  This is automatically determined from `phi_names`.
#' - `trans`/`itrans`: `psi = trans(phi)` and `phi = itrans(psi)`.
#' - `acf(theta, dt, N)`: The increment autocorrelation function in the original basis $\tth$.
#' - `drift(theta, dt, N)`: An optional function specifying the drift on the increment scale.  Thus, linear drift for the trajectory is constant drift for the increments.
#'
#' ## Drift Specification
#'
#' As long as the drift is one of the preset values `none`, `linear`, or `quadratic`, it can be easily specified by the constructor and changed after instantiation, because these methods automatically determine the number of drift parameters `p`.  If a custom drift is used, then not only does it need to specify `p`, but also any additional drift parameters must be added to the acf function as well, such that `n_phi`, `phi_names`, etc. must also be changed.
#'
#' So, it seems that adding a custom drift after instantiation potentially requires changing many members and methods.  So it's probably best not to provide this functionality to users.  If they have a custom drift in mind, they can supply it to the constructor.
#'
#' In this sense, perhaps the base class constructor can have a drift and a `p` argument.  The latter is ignored when drift is one of the default options, but required when drift is a custom function.  Or, the transformed parameter `psi = rep(0, n_phi)` is expected to be valid, in which case a drift evaluation will determine `p`.
#'
#' ## Other Methods

#+ methods
# basics
mod$acf(phi) # original parametrization
mod$drift(phi)
mod$mu_hat(phi)
mod$Sigma_hat(phi)
# change data
mod$dt <- dt
mod$dX <- dX # reallocate memory only if needed

# parameter conversions
mod$trans(phi)
mod$itrans(psi)
mod$trans_full(phi, mu, Sigma)
mod$itrans_full(eta) # returns list

# model fitting
mod$nlp(psi) # computational parametrization
mod$nlp_grad(psi) # if applicable
mod$fisher(eta) # observed fisher information
mod$fit(psi0, vcov = TRUE) # default fitting function???

# model residuals
mod$resid(phi, mu, Sigma) # wrapper to csi_resid()

# simulation
mod$sim(phi, mu, Sigma) # wrapper to csi_sim()

# for AIC and otherwise potentially useful
mod$loglik(phi, mu, Sigma) # original scale

#' # Base Class

#+ base_class
require(R6)

csi_model <- R6Class(
  #' Base class for CSI models.
  #'
  #' @details See helper functions below.
  classname = "csi_model",

  private = list(
    dX_ = NULL, # internal increments.
    dt_ = NULL, # internal interobservation time.
    N_ = NULL, # internal number of increments.
    n_dims = NULL, # internal number of dimensions.
    n_drift = NULL, # internal number of drift basis coefficients.
    n_phi = NULL, # internal number of drift + acf parameters.

    #' Toeplitz matrix object.
    Tz_ = NULL,

    #' Calculate the profile likelihood sufficient statistics.
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
    #' Trajectory increments.
    #'
    #' @details Only reallocates memory for the internal Toeplitz matrix if necessary.
    dX = function(value) {
      if(missing(value)) {
        return(private$dX_)
      } else {
        if(!is.numeric(value) || !is.matrix(value)) {
          stop("`dX` must be a numeric matrix.")
        }
        if(nrow(value) != private$N_) {
          # reallocate Toeplitz matrix
          private$N_ <- nrow(value)
          priviate$n_dims <- ncol(value)
          private$Tz_ <- Toeplitz$new(N = private$N_)
        }
        private$dX_ <- value
      }
    },

    #' Interobservation time.
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
    #' Methods and members to be defined by the derived class.
    trans = NULL,
    itrans = NULL,
    acf = NULL,
    drift = NULL,
    phi_names = NA, # setting to NULL means n_phi = 0

    #' Convert from original to computational basis.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param mu Drift coefficients.
    #' @param Sigma Scale matrix.
    #' @return Full parameter vector in the computational basis.
    trans_full = function(phi, mu, Sigma) {
      c(self$trans(phi), mu, trans_Sigma(Sigma))
    },

    #' Convert from computational to original basis.
    #'
    #' @param eta Vector of parameters in the computational basis.
    #' @return List with elements `phi`, `mu`, and `Sigma`.
    itrans_full = function(eta) {
      n_phi <- private$n_phi
      n_drift <- private$n_drift
      n_dims <- private$n_dims
      # number of parameters in the log-cholesky factor
      n_chol <- n_dims * (n_dims+1) / 2
      phi <- self$itrans(eta[1:n_phi])
      mu <- matrix(eta[n_phi + 1:(n_drift*n_dims)], n_drift, n_dims)
      Sigma <- itrans_Sigma(eta[n_phi + n_drift*n_dims + 1:n_chol])
      list(phi = phi, mu = mu, Sigma = Sigma)
    },

    #' Negative profile likelihood.
    #'
    #' @param psi Kernel parameters in the computational basis.
    nlp = function(psi) {
      phi <- self$itrans(psi) # convert psi to original scale
      suff <- private$get_suff(phi) # sufficient statistics
      -lmn_prof(suff)
    },

    #' Conditional MLE of nuisance parameters.
    #'
    #' @param phi Kernel parameters in the original basis.
    nu_hat = function(phi) {
      suff <- private$get_suff(phi) # sufficient statistics
      list(mu = suff$Bhat, Sigma = suff$S/suff$n)
    },

    #' Loglikelihood function.
    #'
    #' @param phi Kernel parameters in the original basis.
    #' @param mu Drift coefficients.
    #' @param Sigma Scale matrix.
    loglik = function(phi, mu, Sigma) {
      suff <- private$get_suff(phi) # sufficient statistics
      lmn_loglik(Beta = mu, Sigma = Sigma, suff = suff)
    },

    #' Calculate model residuals.
    #'
    #'

    #' Observed Fisher information matrix.
    #'
    #' @param eta Full set of parameters in the computational basis.
    fisher = function(eta) {
      numDeriv::hessian(x = eta, func = function(eta) {
        # convert eta to phi, mu, Sigma
        theta <- self$itrans_full(eta)
        -self$loglik(theta$phi, theta$mu, theta$Sigma)
      })
    },

    #' Model object constructor.
    #'
    #' @param dX Trajectory increments.
    #' @param dt Interobservation time.
    #' @param drift Drift specification.  Either one of the strings "none", "linear", "quadratric", or a function with signature `function(phi, dt, N)`.
    #' @param n_drift Integer number of drift terms.  Ignored if `drift` is one of the default strings.  Required otherwise.
    #'
    #' @details
    #' - Is it worth checking whether model object is valid at construction time?
    #' - `n_phi` is automatically determined from `phi_names`.  But this means the latter must be set by `derived$initialize()` before `super$initialize()` is called.  Should an error be raised if it is not?
    initialize = function(dX, dt, drift = "linear", n_drift) {
      if(!missing(dX)) self$dX <- dX
      self$dt <- dt
      # get drift
      if(drift == "linear") {
        self$drift <- drift_linear
        private$n_drift <- 1
      } else if(self$drift == "none") {
        self$drift <- drift_none
        private$n_drift <- 0
      } else if(self$drift == "quadratic") {
        self$drift <- drift_quadratic
        private$n_drift <- 2
      } else {
        check_drift(drift)
        self$drift <- drift
        private$n_drift <- n_drift
      }
      if(is.na(self$phi_names)) stop("`self$phi_names` has not been set.")
      private$n_phi <- length(self$phi_names)
    }
  )
)

#' Helper functions.
#'
drift_linear <- function(phi, dt, N) dt
drift_none <- function(phi, dt, N) 0
drift_quadratic <- function(phi, dt, N) cbind(dt, diff((0:N*dt)^2))
check_drift <- function(drift) {
  # FIXME: is this too strict?
  if(!identical(methods::formalArgs(drift), c("phi", "dt", "N"))) {
    stop("drift must have argument signature phi, dt, N.")
  }
}


#' Derived class for the FBM model.
fbm_model <- R6Class(
  classname = "fbm_model",
  inherit = csi_model,

  public = list(
    phi_names = "alpha",
    acf = function(phi, dt, N) fbm_acf(alpha = phi, dt = dt, N = N),
    trans = function(phi) logit(phi, 0, 2),
    itrans = function(psi) ilogit(psi, 0, 2)
  )

)

#' Derived class for the fARMA model.
farma_model <- R6Class(
  classname = "farma_model",
  inherit = "csi_model",

  private = list(
    p_ = NULL, # number of AR terms.
    q_ = NULL, # number of MA terms.
    m_ = NULL, # number of terms in MA approximation.

    #' Obtain the `p` and `q` parameters.
    #'
    #' Gives an error if these are not formatted correctly.
    get_pq = function(p, q) {
      if(missing(p)) p <- 0
      if((p - as.integer(p) != 0) && p < 0) {
        "AR order term `p` must be a nonnegative integer."
      }
      if(missing(q)) q <- 0
      if((q - as.integer(q) != 0) && q < 0) {
        "MA order term `q` must be a nonnegative integer."
      }
      setNames(c(p, q), nm = c("p", "q"))
    }
  ),

  public = list(

    trans = function(phi) {
      psi <- phi
      psi[1] <- logit(phi[1], min = 0, max = 2)
      psi[-1] <- logit(phi[-1], min = -1, max = 1)
      psi
    },
    itrans = function(psi) {
      phi <- psi
      phi[1] <- ilogit(psi[1], min = 0, max = 2)
      phi[-1] <- ilogit(psi[-1], min = -1, max = 1)
      phi
    },

    initialize = function(dX, dt, p, q, m = 30, drift = "linear", n_drift) {
      pq <- private$get_pq(p, q)
      p <- pq["p"]
      q <- pq["q"]
      private$p_ <- p
      private$q_ <- q
      private$m_ <- m
      phi_names <- "alpha"
      if(p > 0) phi_names <- c(phi_names, paste0("phi", 1:p))
      if(q > 0) phi_names <- c(phi_names, paste0("rho", 1:q))
      self$phi_names <- phi_names
      self$acf <- function(phi, dt, N) {
        alpha <- phi[1]
        ar_coef <- if(private$p_ == 0) numeric() else phi[1+1:p]
        ma_coef <- if(private$q_ == 0) numeric() else phi[1+p+1:q]
        farma_acf(alpha = alpha, phi = ar_coef, rho = ma_coef,
                  dt = dt, N = N, m = private$m_)
      }
      super$initialize(dX = dX, dt = dt, drift = drift, n_drift = n_drift)
      ## FIXME: modify the drift term
    }
  )
)

fobj <- fbm_model$new()

fobj$foo(5)

#' Create a CSI derived class generator.
csi_model <- function(classname,
                      trans,
                      itrans,
                      acf,
                      drift = "linear") {
  if(drift == "linear") {
    drift <- function(phi, dt, N) dt
  } else if(drift == "none") {
    drift <- function(phi, dt, N) 0
  } else if(drift == "quadratic") {
    drift <- function(phi, dt, N) cbind(dt, diff((0:N*dt)^2))
  }
  if(missing(trans) != missing(itrans)) {
    stop("`trans` and `itrans` must both be missing or both be provided.")
  }
  if(missing(trans)) trans <- function(phi) phi
  if(missing(itrans)) itrans <- function(psi) psi
  R6Class(
    classname = classname,
    inherit = csi_model,
    public
  )
}

#--- tests ---------------------------------------------------------------------

require(R6)

base_class <- R6Class(
  "base_class",
  private = list(
    x_= 10
  ),
  ## active = list(
  ##   x = function() private$x_
  ## )
)

derived_class <- R6Class(
  "derived_class",
  inherit = base_class,
  private = list(
    x_ = 20
  ),
  ## active = list(
  ##   x = function() private$dx_
  ## ),
  public = list(
    get_x = function(base = FALSE) {
      if(base) return(super$private$x_)
      private$x_
    }
  )
)

foo <- derived_class$new()

foo$get_x()
foo$get_x(base=TRUE) # fails
