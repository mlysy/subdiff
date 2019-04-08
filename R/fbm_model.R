#' Constructor of fBM model object.
#' 
#' @details Constructor function that generates a list of fBM model for purpose of fitting.
#' @template ret-csi_class
#' 
#' @example examples/fbm_model.R
#' 
#' @export
fbm_model <- function() {
  fbm_theta <- "alpha"
  fbm_trans <- function(alpha) logit(alpha, min = 0, max = 2)
  fbm_itrans <- function(theta) ilogit(theta, min = 0, max = 2)
  fbm_penalty <- function(gamma) 0
  fbm_penalty <- function(gamma) 0
  
  model <- list(theta_names = fbm_theta,
                acf = fbm_acf,
                theta_trans = fbm_trans,
                theta_itrans = fbm_itrans,
                penalty = fbm_penalty)
  class(model) <- "csi_model"  
  model
}
