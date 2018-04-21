#' @title Downsample Function
#' @param Yt Time series
#' @param ds Downsampling rate
#' @param pos Position index
#' @export
downSample <- function(Yt, ds, pos = 1) {
  if(pos > ds) {
    stop("position index must be smaller than downsample rate.")
  }
  if(ds == 1) {
    Yt
  } else {
    N <- nrow(Yt)
    N.ds <- floor(N/ds)
    Yt.ds <- Yt[seq(from = pos, by = ds, length.out = N.ds), ]
    Yt.ds
  }
}
