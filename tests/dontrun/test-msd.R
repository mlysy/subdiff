#--- test the numerical msd ----------------------------------------------------

require(subdiff)
require(testthat)

# fbm model

fbm_obj <- fbm_model$new(dX = matrix(0), dt = 1)

n <- sample(100, 1)
alpha <- runif(1, 0, 2)
dt <- runif(1)
t <- runif(n, -100, 100)
## t <- dt * 0:n

msd1 <- fbm_obj$msd(phi = alpha, t = t[3])
msd2 <- fbm_msd(t = t[3], alpha = alpha)
plot(t, msd1, type = "l", log = "xy")
lines(t, msd2, col = "red")
expect_equal(log(msd1), log(msd2), tolerance = 1e-3)

#--- test msd vs acf -----------------------------------------------------------

# fbm model
alpha <- runif(1, 0, 2)
t <- runif(1, -100, 100)
expect_equal(fbm_msd(t, alpha), fbm_acf(alpha, abs(t), 1))

# fsd model
dt <- runif(1)
alpha <- runif(1, 0, 2)
tau <- runif(1)
sigma2 <- runif(1)
t <- runif(1)
# in acf mode, dt = abs(t), so need to rescale tau
tau2 <- tau * dt / abs(t)
expect_equal(fsd_msd(t, dt, alpha, tau, sigma2),
             fsd_acf(alpha, tau2, sigma2, abs(t), 1))

#--- ok check fsd integral formula ---------------------------------------------

require(cubature)

# covariance of fbm
fbm_cov <- function(t, s, alpha) {
  .5 * (abs(t)^alpha + abs(s)^alpha - abs(t-s)^alpha)
}

gfun <- function(t, tau, alpha) {
  out <- abs(t+tau)^(alpha+2) + abs(t-tau)^(alpha+2) - 2*abs(t)^(alpha+2)
  out/(2*tau^2*(alpha+1)*(alpha+2))
}

hfun <- function(t, tau, alpha) {
  -(abs(t-tau)^(alpha+1) - t^(alpha+1)) / (2*tau*(alpha+1))
}

fsd_cov <- function(t, s, tau, alpha) {
  hfun(t, tau, alpha) + hfun(s, tau, alpha) - gfun(abs(t-s), tau, alpha)
}

alpha <- runif(1)
t <- runif(1)
s <- runif(1)
tau <- runif(1) * min(t, s) # absolute units


pcubature(
  f = function(x) fbm_cov(t-x[1], s-x[2], alpha)/tau^2,
  ## f = function(x) abs(t - x[1])^alpha/(2*tau^2),
  lowerLimit = c(0, 0),
  upperLimit = c(tau, tau),
  fDim = 1
)

## hfun(t, tau, alpha)

fsd_cov(t, s, tau, alpha)


# ok now check fsd msd again

# fsd model
dt <- runif(1)
alpha <- runif(1, 0, 2)
tau <- runif(1)
sigma2 <- runif(1)
tau2 <- tau*dt # absolute units

t <- dt + runif(1) # offset
fsd_cov(t+dt, t+dt, tau2, alpha) +
  fsd_cov(t, t, tau2, alpha) -
  2 * fsd_cov(t, t+dt, tau2, alpha)
2 * (gfun(dt, tau2, alpha) - gfun(0, tau2, alpha))


fsd_acf(alpha, tau, 0, dt, 1)
fsd_msd(dt, dt, alpha, tau, 0)

t <- 2 * dt # runif(1, -100, 100)

# in acf mode, dt = abs(t), so need to rescale tau
tau2 <- tau * dt / abs(t)
fsd_acf(alpha, tau2, sigma2, abs(t), 1)
2*hfun(t, tau2, alpha) - gfun(0, tau2, alpha)

expect_equal(fsd_msd(t, dt, alpha, tau, sigma2),
             fsd_acf(alpha, tau2, sigma2, abs(t), 1))
