
# dimension of covariacne matrix
getq <- function(ndims) if(ndims == 1) 1 else 3

# max of min of abs and rel error
max.xdiff <- function(x) {
  xdiff <- abs(diff(x))
  max(pmin(xdiff[,1], xdiff[,2]))
}

# log(1+exp(x))
log1pexp <- function(x) {
  n <- length(x)
  y <- rep(NA, n)
  ind <- x <= -37
  if(any(ind)) y[ind] <- exp(x[ind])
  ind <- x > -37 & x <= 18
  if(any(ind)) y[ind] <- log1p(exp(x[ind]))
  ind <- x > 18 & x <= 33.3
  if(any(ind)) y[ind] <- x[ind] + exp(-x[ind])
  ind <- x > 33.3
  if(any(ind)) y[ind] <- x[ind]
  y
}
