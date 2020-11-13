#' Constructor of fsd model object.
#'
#' @details Constructor function that generates a list of fsd model for purpose of fitting.
#' @template ret-csi_class
#'
#' @example examples/fsd_model.R
#'
#' @export
fsd_model <- function() {
  # parameter name
  fsd_theta <- c("alpha", "tau", "sigma2")

  # ACF function
  .fsd_acf <- function(theta, dt, N) {
    fsd_acf(theta[1], theta[2], theta[3], dt, N)
  }

  # transformation, from original scale (theta) to unrestricted scale (gamma)
  fsd_trans <- function(theta) {
    gamma <- theta
    # alpha
    gamma[1] <- logit(theta[1], min = 0, max = 2)
    # tau
    gamma[2] <- logit(theta[2], min = 0, max = 1)
    # sigma2
    gamma[3] <- log(theta[3])
    gamma
  }

  # inverse transformation, from unrestricted scale (gamma) to original scale (theta)
  fsd_itrans <- function(gamma) {
    theta <- gamma
    # alpha
    theta[1] <- ilogit(gamma[1], min = 0, max = 2)
    # tau
    theta[2] <- ilogit(gamma[2], min = 0, max = 1)
    # sigma2
    theta[3] <- exp(gamma[3])
    theta
  }

  # penalty function on transformed scale
  fsd_penalty <- function(gamma) {
    log1pe(gamma[2]) + log1pe(-gamma[2]) - gamma[3]
  }

  # plug-in values
  # fsd_plugin <- names_plugin <- NULL
  # if(!missing(tau)) {
  #   fsd_plugin <- c(fsd_plugin, tau)
  #   names_plugin <- c(names_plugin, "tau")
  # }
  # if(!missing(sigma2)) {
  #   fsd_plugin <- c(fsd_plugin, sigma2)
  #   names_plugin <- c(names_plugin, "sigma2")
  # }
  # names(fsd_plugin) <- names_plugin

  model <- list(theta_names = fsd_theta,
                acf = .fsd_acf,
                theta_trans = fsd_trans,
                theta_itrans = fsd_itrans,
                penalty = fsd_penalty)
  class(model) <- "csi_model"
  model
}
