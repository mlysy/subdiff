#' Constructor of fLOC model object.
#' 
#' @details Constructor function that generates the default S3 class of fLOC model.
#' @template ret-csi_class
#' 
#' @examples 
#' model <- floc_model()
#' 
#' @export
floc_model <- function() {
  # parameter name
  floc_theta <- c("alpha", "tau", "sigma2")
  
  # ACF function
  .floc_acf <- function(theta, dT, N) {
    floc_acf(theta[1], theta[2], theta[3], dT, N)
  }
  
  # transformation, from original scale (theta) to unrestricted scale (gamma)
  floc_trans <- function(theta) {
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
  floc_itrans <- function(gamma) {
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
  floc_penalty <- function(gamma) {
    log1pe(gamma[2]) + log1pe(-gamma[2]) - gamma[3]
  }
  
  # plug-in values
  # floc_plugin <- names_plugin <- NULL
  # if(!missing(tau)) {
  #   floc_plugin <- c(floc_plugin, tau)
  #   names_plugin <- c(names_plugin, "tau")
  # }
  # if(!missing(sigma2)) {
  #   floc_plugin <- c(floc_plugin, sigma2)
  #   names_plugin <- c(names_plugin, "sigma2")
  # }
  # names(floc_plugin) <- names_plugin
  
  model <- list(theta_names = floc_theta,
                acf = .floc_acf,
                theta_trans = floc_trans,
                theta_itrans = floc_itrans,
                penalty = floc_penalty)
  class(model) <- "csi_model"  
  model
}
