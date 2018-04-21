# multivariate cochran test

k <- 10
p <- 3
rMnorm <- function(n,p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

est <- rMnorm(k,p)
ve <- replicate(k, crossprod(rMnorm(p)))

# calculate by hand
est0 <- 0
ve0 <- 0
for(ii in 1:k) {
  est0 <- est0 + solve(ve[,,ii], est[ii,])
  ve0 <- ve0 + solve(ve[,,ii])
}
ve0 <- solve(ve0)
est0 <- ve0 %*% est0
QQ <- 0
for(ii in 1:k) {
  QQ[1] <- QQ + crossprod(est[ii,]-est0, solve(ve[,,ii], est[ii,]-est0))
}

Qt <- cochranMQ(est, ve)

est0-Qt$est
ve0-Qt$ve
QQ-Qt$Q

# ok check the null distribution for fun

n <- 10000
k <- 10
p <- 3

nreps <- 1000
Qt <- replicate(nreps, {
  X <- lapply(1:k, function(ii) rMnorm(n,p))
  est <- as.matrix(sapply(X, colMeans))
  if(p > 1) est <- t(est)
  ve <- array(sapply(X, var), dim = c(p,p,k))/n
  cochranMQ(est = est, ve = ve)$Q
})

Qt <- replicate(nreps, {
  X <- rMnorm(n,k)
  xbar <- colMeans(X)
  S <- apply(X, 2, var)
  cochranQ(xbar, sqrt(S/n))["pval"]
})

hist(Qt, breaks = 50, freq = FALSE)
curve(dchisq(x, df = p*(k-1)), add = TRUE, col = "red")
