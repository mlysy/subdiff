#' Normalizing transformation for the AR/MA filters.
#'
#' @param rho Filter parameter on regular scale.
#' @return Filter parameter on the regular or normalized scale.
#' @details The normalizing transformation for the filter parameter is
#' \deqn{
#' \eta = logit((\rho+1)/2).
#' }
#' @name trans_rho
#' @export
trans_rho <- function(rho) log(1+rho) - log(1-rho)

#' @rdname trans_rho
#' @param eta Filter parameter on the normalized scale.
#' @export
itrans_rho <- function(eta) 2/(1+exp(-eta)) - 1
