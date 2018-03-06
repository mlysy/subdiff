#' Normalizing transformation for the dynamic error model.
#'
#' @param tau Camera exposure time divided by the interobservation time.
#' @return Camera exposure parameter on the regular or normalized scale.
#' @details The normalizing transformation for the camera exposure parameter is
#' \deqn{
#' \tau = logit(\tau).
#' }
#' @name trans_tau
#' @export
trans_tau <- function(tau) log(tau/(1-tau))

#' @rdname trans_tau
#' @param tau Camera exposure parameter on the normalized scale.
#' @export
itrans_tau <- function(tau) 1/(1+exp(-tau))
