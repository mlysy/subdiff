#' Calculate the mean square displacement of the fBM model.
#'
#' @template args-t
#' @param alpha Subdiffusion exponent. A scalar between 0 and 2.
#'
#' @template ret-msd
#'
#' @details If `X_t` is an fBM process, then
#' ```
#' msd_X(t) = t^alpha.
#' ```
#' @export
fbm_msd <- function(t, alpha) abs(t)^alpha
