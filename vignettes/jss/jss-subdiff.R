#--- R scripts for JSS paper ---------------------------------------------------

# load required packages
require(subdiff)
require(mniw)
require(LMN)
require(SuperGauss)
require(parallel)
require(numDeriv)
require(ggplot2)
require(kableExtra)

#--- Load HBE dataset and plot trajectories ------------------------------------

# FIXME: include all mucus concentrations?
head(hbe) # data is already part of subdiff package
npaths <- length(unique(hbe$id)) # number of trajectories
n_rng <- range(tapply(hbe$id, hbe$id, length)) # range of trajectory lengths
dT <- 1/60 # interobservation time

# plot trajectories and histogram of trajectory lengths
qplot(x = x, y = y, data = hbe,
      color = factor(id), geom = "line",
      xlab = expression("X Position ("*mu*m*")"),
      ylab = expression("Y Position ("*mu*m*")")) +
  theme(legend.position = "none")
qplot(x = tapply(hbe$id, hbe$id, length), geom = "histogram",
      xlab = "Number of Observations", ylab = "Counts")


#--- Empirical MSDs ------------------------------------------------------------

# extract variables from a dataset
get_vars <- function(data, id, vars) data[data$id == id, vars]

# calculate empirical MSDs
ids <- unique(hbe$id)
nlags <- 500 # number of lags
dT <- 1/60 # interobservation time
tseq <- 1:nlags*dT # time sequence
msd_hat <- sapply(ids, function(id) {
  Xt <- get_vars(data = hbe, id = id, vars = c("x", "y"))
  cbind(id = id, t = tseq, msd = msd_fit(Xt = Xt, nlag = nlags))
}, simplify = FALSE)
msd_hat <- as.data.frame(do.call(rbind, msd_hat))

# calculate mean and standard error
msd_stats <- tapply(msd_hat$msd, msd_hat$t, function(msd) {
  c(est = mean(msd), se = sd(msd)/sqrt(length(ids)))
})
names(msd_stats) <- NULL
msd_stats <- as.data.frame(cbind(t = tseq, do.call(rbind, msd_stats)))

# plot msd's with mean +/- 1.96 se
ggplot(data = msd_stats) + # blank plot
  # individual msd's
  geom_line(data = msd_hat,
            mapping = aes(x = t, y = msd, group = factor(id),
                          color = "msd")) +
  # +/- 1.96*se confidence band
  geom_ribbon(mapping = aes(x = t,
                            ymin = est-1.96*se, ymax = est+1.96*se,
                            color = "se", fill = "se", alpha = "se")) +
  # mean line
  geom_line(aes(x = t, y = est, color = "est")) +
  # set colors manually
  scale_color_manual(values = c(msd = "lightblue", est = "black", se = NA)) +
  scale_fill_manual(values = c(msd = NA, est = NA, se = "red")) +
  scale_alpha_manual(values = c(msd = NA, est = NA, se = .5)) +
  # axes names and log-scale
  scale_x_continuous(name = expression("Time (s)"), trans = "log10") +
  scale_y_continuous(name = expression("MSD ("*mu*m^2*")"), trans = "log10") +
  # remove legend
  theme(legend.position = "none")

#--- estimate & plot subdiffusion parameters -----------------------------------

aD_hat <- sapply(ids, function(id) {
  msd <- get_vars(data = msd_hat, id = id, vars = "msd")
  aD <- msd_ls(msd = msd, tseq = tseq)
  c(id = id, aD)
}, simplify = FALSE)
aD_hat <- as.data.frame(do.call(rbind, aD_hat))

qplot(x = D, y = alpha, data = aD_hat, group = factor(id), log = "x",
      xlab = expression(hat(D)), ylab = expression(hat(alpha)))


#--- simulate fBM trajectories -------------------------------------------------

# true parameter values
alpha <- .8 # subdiffusion exponent
mu <- c(.025, -.025) # drift coefficients
Sigma <- cbind(c(.1,.05), c(.05,.1)) # scale matrix

N <- 1800 # number of observations
dT <- 1/60 # interobservation time

# simulate an fBM trajectory

# step 1: obtain fbm covariance matrix
fbm_cov <- function(t, s, alpha) {
  .5 * (abs(t)^alpha + abs(s)^alpha - abs(t-s)^alpha)
}
tseq <- (1:N)*dT # sequence of time points excluding t = 0
V <- outer(tseq, tseq, fbm_cov, alpha = alpha) # temporal variance matrix

# step 2: simulate fbm
system.time({
  Xt <- matrix(rnorm(N*2), N, 2) # generate white noise
  Xt <- t(chol(V)) %*% Xt # temporal correlations
  Xt <- Xt %*% chol(Sigma) # spatial correlations
  Xt <- Xt + tseq %o% mu # drift
  Xt <- cbind(0, Xt) # add first observation
})

# simulate fbm trajectory efficiently
system.time({
  acf <- fbm_acf(alpha, dT = dT, N = N) # increment autocorrelation
  # generate temporally correlated increments
  dX <- SuperGauss::rSnorm(n = 2, acf = acf)
  dX <- dX %*% chol(Sigma) # spatial correlations
  dX <- dX + rep(dT, N) %o% mu # drift
  Xt <- apply(dX, 2, cumsum) # convert increments to positions
  Xt <- cbind(0, Xt) # add first observations
})


#--- fit fBM to all trajectories in a given data frame -------------------------

ids <- unique(hbe$id) # trajectory labels
dT <- 1/60 # interobservation time

# create a function to fit the MLE for dataset id
# FIXME: should wrap fbm_fit in a try-catch block
# should also save to disk after each fit
fbm_mle <- function(id) {
  # extract the dataset
  Xt <- get_vars(data = hbe, id = id, vars = c("x", "y"))
  # convert to increments
  dX <- apply(Xt, 2, diff)
  N <- nrow(dX) # number of increments
  # external memory allocation (to speed up calculations)
  Tz <- SuperGauss::Toeplitz(n = N)
  fbm_fit(dX = dX, dT = dT, Tz = Tz, var_calc = TRUE)
}

parallel::detectCores() # detect number of cores
# on my system half the cores are virtual and don't increase performance much
ncores <- 4

# create parallel cluster
parclust <- parallel::makeCluster(ncores)
# add required package(s) to cluster
parallel::clusterEvalQ(cl = parclust, expr = {
  require(subdiff)
})
# add required global variables to cluster
parallel::clusterExport(cl = parclust,
                        varlist = c("get_vars", "hbe", "dT"))

# fit MLE to each dataset in parallel
system.time({
  Psi_est <- parSapply(cl = parclust, X = ids, FUN = fbm_mle)
})

# shut down cluster to free up resources
parallel::stopCluster(cl = parclust)


# boxplots for MLE of each parameter
Psi_mle <- as.data.frame(cbind(id = ids, do.call(rbind, Psi_est["coef",])))
# rearrange to have parameter name and value as a variables
param_names <- c("kappa", "mu[1]", "mu[2]",
                 "lambda[1]", "lambda[2]", "lambda[3]")
names(param_names) <- names(Psi_mle)[-1]
Psi_mle <- do.call(rbind, lapply(names(param_names), function(par) {
  data.frame(id = ids, parameter = as.character(param_names[par]),
             value = Psi_mle[,par])
}))
# plot
qplot(y = value, data = Psi_mle, # empty plot
      geom = "blank", xlab = "", ylab = "", xlim = .5*c(-1,1)) +
  # horizontal whisker bars in background
  stat_boxplot(geom = "errorbar", width = .25) +
  # boxplot in foreground
  geom_boxplot() +
  # separate into different plots with different axes
  facet_wrap(facets = ~ parameter, scales = "free",
             labeller = label_parsed) +
  # remove x axis
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


#--- confidence intervals ------------------------------------------------------

# contents of model_fit object
# FIXME: write a `show` method?
# FIXME: standardize the computational basis
sapply(Psi_est[,1], signif, digits = 3)

# Get log-nondimensional D (log_DN) from (alpha, logD),
# the bead diameter (d) and the diffusivity of the bead in water (Dw).
nondim <- function(alpha, logD, d = 1, logDw = log(.46)) {
  (2 * (alpha-1)) * log(d) + alpha * (logD - logDw)
  ## (d^(2*(alpha-1)))*(D/(Dw)^alpha)
}

# transform psi = (kappa, lambda1) to
# omega = (alpha, D, log10 D, D_N, log10 D_N)
Gfun <- function(psi) {
  log10e <- log10(exp(1)) # unit conversion between log_e and log_10
  alpha <- ilogit(psi[1], 0, 2)
  logD <- psi[2] - log(4)
  logDN <- nondim(alpha, logD)
  onames <- c("alpha", "D", "logD", "DN", "logDN")
  setNames(c(alpha,
             exp(logD), log10e * logD,
             exp(logDN), log10e * logDN), nm = onames)
}

# FIXME: change `aD` to `omega`, and possibly a few more issues now that
# we have 5 outputs instead of 2.
# for psi = (kappa, lambda1) and omega = (alpha, logD), calculate:
# - psi_mle, psi_se
# - aD_mle, aD_se
# - aD_cov = cov(alpha, logD)
aD_stat <- sapply(ids, function(id) {
  # MLE and variance estimator on unconstrained scale (kappa, lambda1)
  ipsi <- c(1,4) # indices of (kappa, lambda1) in psi
  psi_mle <- Psi_est[["coef",id]][ipsi]
  psi_ve <- Psi_est[["vcov",id]][ipsi, ipsi]
  # MLE of aD
  aD_mle <- Gfun(psi_mle)
  # variance estimate of aD_mle
  # jacobian
  # NOTE: this is the transpose of how it is defined in the vignette
  Jac <- numDeriv::jacobian(func = Gfun, x = psi_mle)
  aD_ve <- tcrossprod(Jac %*% psi_ve, Jac) # t(Jac) %*% psi_ve %*% Jac
  # standard errors of aD_mle
  aD_se <- sqrt(diag(aD_ve))
  # cov(alpha_mle, logD_mle)
  aD_cov <- aD_ve[1,2]
  setNames(c(id, psi_mle, sqrt(diag(psi_ve)), aD_mle, aD_se, aD_cov),
           c("id", "kappa_mle", "lambda1_mle",
             "kappa_se", "lambda1_se",
             "alpha_mle", "logD_mle",
             "alpha_se", "logD_se", "aD_cov"))
}, simplify = FALSE)
aD_stat <- as.data.frame(do.call(rbind, aD_stat))

# calculate confidence intervals for alpha using two methods:

# method 1: normal approximation on alpha_mle
alpha_ci1 <- cbind(L = aD_stat$alpha_mle - 1.96 * aD_stat$alpha_se,
                   U = aD_stat$alpha_mle + 1.96 * aD_stat$alpha_se)
# method 2: normal approximation on kappa_mle
alpha_ci2 <- cbind(L = ilogit(aD_stat$kappa_mle - 1.96 * aD_stat$kappa_se, 0, 2),
                   U = ilogit(aD_stat$kappa_mle + 1.96 * aD_stat$kappa_se, 0, 2))

# compare 2 methods via table
tab <- cbind(aD_stat[,"id"], alpha_ci1, alpha_ci2)
tab <- signif(tab[1:10,], digits = 3)
colnames(tab) <- c("Trajectory ID", "Lower", "Upper", "Lower", "Upper")
require(kableExtra) # for nice HTML tables
kable(tab) %>%
  kable_styling(bootstrap_options = c("striped", "responsive"), full_width = TRUE) %>%
  add_header_above(c(" " = 1, "Direct Method" = 2, "Delta Method" = 2)) %>%
  column_spec(1:5, width = "3cm")


#--- display confidence ellipses -----------------------------------------------

# points of a 2D ellipse
ellipse <- function (mu, V, alpha = 0.95, n = 100) {
  eV <- eigen(V)
  hlen <- sqrt(qchisq(alpha, df = 2) * eV$val)
  phi <- atan2(eV$vec[2, 1], eV$vec[1, 1])
  theta <- seq(0, 2 * pi, len = n + 1)
  x <- hlen[1] * cos(theta)
  y <- hlen[2] * sin(theta)
  alpha <- atan2(y, x)
  rad <- sqrt(x^2 + y^2)
  cbind(x = rad * cos(alpha + phi) + mu[1],
        y = rad * sin(alpha + phi) + mu[2])
}

# calculate confidence ellipses
aD_ell <- sapply(ids, function(id) {
  # extract aD_mle and aD_ve
  aD_mle <- as.numeric(get_vars(aD_stat, id = id,
                                vars = c("alpha_mle", "logD_mle")))
  aD_ve <- matrix(NA, 2, 2)
  diag(aD_ve) <- as.numeric(get_vars(aD_stat, id = id,
                                     vars = c("alpha_se", "logD_se")))^2
  aD_ve[1,2] <- aD_ve[2,1] <- get_vars(aD_stat, id = id, vars = "aD_cov")
  # calculate points of 95% confidence ellipse
  ell <- ellipse(mu = aD_mle, V = aD_ve)
  colnames(ell) <- c("alpha", "logD")
  cbind(id = id, ell)
}, simplify = FALSE)
aD_ell <- as.data.frame(do.call(rbind, aD_ell))

# plot
ggplot(data = aD_ell) + # empty plot
  # confidence ellipses
  geom_polygon(aes(x = logD, y = alpha, group = id,
                   fill = "ci", alpha = "ci", color = "ci")) +
  # mles
  geom_point(data = aD_stat,
             mapping = aes(x = logD_mle, y = alpha_mle, group = id,
                           color = "mle")) +
  # custom colors
  scale_fill_manual(values = c(mle = "black", ci = "red")) +
  scale_color_manual(values = c(mle = "black", ci = NA)) +
  scale_alpha_manual(values = c(mle = 0, ci = .2)) +
  # axis labels
  scale_x_continuous(name = expression(log[10]*(D))) +
  scale_y_continuous(name = expression(alpha)) +
  # remove legend
  theme(legend.position = "none")
