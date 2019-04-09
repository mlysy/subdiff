#' Generalized logistic and inverse-logistic transformations.
#'
#' @name logit
#' @param x A numeric scalar or vector.
#' @param min,max Numeric scalar or vectors denoting the bounds of the transformations.
#' @return The logit or inverse-logit transoformation.
#'
#' @details The generalized logistic transformation is
#' \preformatted{
#' logit(x) = log(p/(1-p)),    p = (x-min)/(max - min),
#' }
#' and its inverse is
#' \preformatted{
#' ilogit(x) = min + (max-min)/(1+exp(-x).
#' }
#' @export
logit <- function(x, min = 0, max = 1) {
  x <- (x - min) / (max - min)
  log(x) - log(1-x)
}

#' @rdname logit
#' @export
ilogit <- function(x, min = 0, max = 1) {
  1/(1+exp(-x)) * (max - min) + min
}
