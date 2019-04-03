#' Normalizing transformation for the subdiffusion coefficient.
#'
#' @param alpha Subdiffusion coefficient on regular scale.
#' @return Subdiffusion coefficient on the regular or normalized scale.
#' @details The normalizing transformation for the subdiffusion coefficient is
#' \deqn{
#' \gamma = logit(\alpha/2).
#' }
#' @name trans_alpha
#' @export
trans_alpha <- function(alpha) log(alpha/(2-alpha))

#' @rdname trans_alpha
#' @param gamma Subdiffusion coefficient on normalized scale.
#' @export
itrans_alpha <- function(gamma) 2/(1+exp(-gamma))

