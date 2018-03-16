
# dimension of covariacne matrix
getq <- function(ndims) if(ndims == 1) 1 else 3

# max of min of abs and rel error
max.xdiff <- function(x) {
  xdiff <- abs(diff(x))
  max(pmin(xdiff[,1], xdiff[,2]))
}
