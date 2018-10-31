#' Normalizing transformation for the subdiffusion coefficient.
#'
#' @param gamma Subdiffusion coefficient on regular scale.
#' @return Subdiffusion coefficient on the regular or normalized scale.
#' @details The normalizing transformation for the subdiffusion coefficient is
#' \deqn{
#' \theta = logit((\gamma+1)/2).
#' }
#' @name trans_gamma
#' @export
trans_gamma <- function(gamma) {
  log((1+gamma) / (1-gamma))
}

#' @rdname trans_gamma
#' @param igamma Subdiffusion coefficient on normalized scale.
#' @export
itrans_gamma <- function(igamma) {
  1 - 2/(exp(igamma) + 1)
}
