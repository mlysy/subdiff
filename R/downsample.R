#' Downsample a multi-dimensional time series.
#'
#' @param Xt One or two-column matrix of observations.
#' @param r Integer of downsampling rate.
#' @param position Integer of position index for downsampled series, smaller than or equal to \code{r}.
#' @return An matrix of downsampled observation.
#' 
#' @details Sample observations from original time series \eqn{X = [X_0, X_2, X_3, ..., X_N]} at constant rate \code{r}. For given position index \code{i} finally returns \eqn{X_{ds} = [X_i, X_{i+r}, X_{i+2r}, ...]}.
#' 
#' @example examples/Xt_setup.R
#' @example examples/downsample.R
#' 
#' @export
downsample <- function(Xt, r, position = 1) {
  if(position > r) stop("position index cannot be larger than downsample rate")
  as.matrix(Xt[seq(from = position, by = r, len = floor(nrow(Xt)/r)), ])
}
