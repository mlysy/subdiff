#' Constructor for farma model object.
#'
#' @param p An integer of Auto-regressive order.
#' @param q An integer of Moving-average order.
#' @details Constructor function that generates a list of farma(p,q) model for purpose of fitting.
#' @template ret-csi_class
#'
#' @example examples/farma_model.R
#'
#' @export
farma_model <- function(p, q) {
  if(!p & !q) {
    # depreciated to fBM model
    return(fbm_model())
  }

  farma_theta <- "alpha"
  if(p) farma_theta <- c(farma_theta, paste0("phi", 1:p))
  if(q) farma_theta <- c(farma_theta, paste0("rho", 1:q))

  .farma_acf <- function(theta, dT, N, m = 30) {
    if(!p) {
      farma_acf(theta[1], numeric(), theta[1+1:q], dT, N, m)
    } else {
      farma_acf(theta[1], theta[1+1:p], theta[1+p+1:q], dT, N, m)
    }
  }

  farma_trans <- function(theta) {
    gamma <- theta
    # alpha
    gamma[1] <- logit(theta[1], min = 0, max = 2)
    # phi and rho
    gamma[-1] <- logit(theta[-1], min = -1, max = 1)
    gamma
  }

  farma_itrans <- function(gamma) {
    theta <- gamma
    # alpha
    theta[1] <- ilogit(gamma[1], min = 0, max = 2)
    # phi and rho
    theta[-1] <- ilogit(gamma[-1], min = -1, max = 1)
    theta
  }

  farma_penalty <- function(gamma) 0

  model <- list(theta_names = farma_theta,
                acf = .farma_acf,
                theta_trans = farma_trans,
                theta_itrans = farma_itrans,
                penalty = farma_penalty)
  class(model) <- "csi_model"
  model
}
