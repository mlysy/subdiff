#' Downsample the observation of time series.
#'
#' @param Xt One or two-column matrix of observations.
#' @param r Integer of downsampling rate.
#' @param position Integer of position index for downsampled series, smaller than or equal to `r`.
#' @return An matrix of downsampled observation.
#' 
#' @details Sample observations from original time series \eqn{X = [X_0, X_2, X_3, ..., X_N]} at constant rate `r`. For given position index `i` finally returns \eqn{X_{ds} = [X_i, X_{i+r}, X_{i+2r}, ...]}.
#' 
#' @examples 
#' Xt <- apply(matrix(rnorm(2000), 2, 1000), 2, cumsum) # generates a 2-d series of Brownian motion.
#' Xt_ds2 <- downsmaple(Xt, r = 2)
#' 
#' @export
downsample <- function(Xt, r, position = 1) {
  if(position > r) stop("position index cannot be larger than downsample rate")
  as.matrix(Xt[seq(from = position, by = r, len = floor(nrow(Xt)/r)), ])
}
