#' @title Downsample Function
#' @param Yt Time series
#' @param ds Downsampling rate
#' @param pos Position index
#' @export

downSample <- function(Yt, ds, pos = 1) {
  if(ds == 1){
    Yt
  } else {
    if(pos > ds) {
      stop("position index must be smaller than downsample rate.")
    }
    N <- nrow(Yt)
    d <- ncol(Yt)
    N.ds <- floor(N / ds)
    Yt.ds <- matrix(NA, N.ds, d)
    for(ii in 1:N.ds) {
      Yt.ds[ii, ] <- Yt[pos + (ii-1) * ds, ]
    }
    Yt.ds
  }
}
