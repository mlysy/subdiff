#--- code for vignettes --------------------------------------------------------

# first, generate some interesting simulated data.
# do this by fitting a few real datasets with the GLE model

require(subdiff)
require(tidyverse)

hbe <- read.csv("/Users/mlysy/Dropbox/Shared/gpsi/HBE data/HBE_data.csv")
hbe <- as_tibble(hbe) %>% filter(wt == "2p5")

# focus on wt = 2p5 as in the paper
hbe %>%
  mutate(x = x + setNames(10*rnorm(n_distinct(localid)),
                          unique(localid))[as.character(localid)],
         y = y + setNames(10*rnorm(n_distinct(localid)),
                          unique(localid))[as.character(localid)]) %>%
  ggplot(aes(x = x, y = y, color = factor(localid))) +
  geom_line()

# ok.  now fit crouse gle to each dataset, and resimulate.

crouse_acf <- function(eta, dT, N) {
  tseq <- dT*1:N
  msd <- crouse_msd(tseq,
                    alpha = eta[1], tmin = eta[2], tmax = eta[3])
  msd2acf(msd = msd)
}

crouse_fit <- function(ii, test_out = FALSE) {
  ## Xt <- hbe %>% filter(localid == 1) %>% select(x, y) %>% as.matrix
  Xt <- as.matrix(hbe[hbe$localid == ii, c("x", "y")])
  dX <- apply(Xt, 2, diff)
  dT <- 1/60
  N <- nrow(dX)
  tseq <- dT*1:N
  ## plot(dT*(1:100), msd_fit(Xt, nlag = 100))
  Tz <- Toeplitz(n = N)
  # sufficient statistics
  crouse_suff <- function(eta) {
    Tz$setAcf(acf = crouse_acf(eta, dT = dT, N = N))
    lmn.suff(Y = dX, X = dT, V = Tz, Vtype = "acf")
  }
  # negative profile likelihood
  crouse_nll <- function(eta) {
    if((abs(eta[1] - .5)) > .5) return(Inf) # 0 < alpha < 1
    if(any(eta[2:3] < 0)) return(Inf) # tmin, tmax > 0
    nlp <- -lmn.prof(suff = crouse_suff(eta))
    # penalty on tmax
    nlp <- nlp - dnorm(log10(eta[3]), mean = log10(30), sd = 3, log = TRUE)
    nlp
  }
  # fit
  eta0 <- c(alpha = .8, tmin = dT, tmax = dT*N)
  opt <- optim(par = eta0, fn = crouse_nll)
  # get mle
  eta_mle <- opt$par
  suff_mle <- crouse_suff(eta_mle)
  mu_mle <- setNames(drop(suff_mle$Bhat), c("mu1", "mu2"))
  Sigma_mle <- suff_mle$S/suff_mle$n
  theta_mle <- c(eta_mle, mu_mle,
                 setNames(Sigma_mle[c(1,2,4)],
                          c("Sigma11", "Sigma21", "Sigma22")))
  # output
  if(test_out) {
    D_mle <- mean(diag(Sigma_mle))
    out <- sapply(ls(), get, envir = -3) # output everything from function
  } else {
    out <- c(id = ii, N = N, theta_mle) # output mle only
  }
  out
}

out <- crouse_fit(1, test_out = TRUE)

Xt <- out$Xt
tseq <- out$tseq
eta_mle <- out$eta_mle
msd_mle <- out$D_mle * crouse_msd(tseq, eta_mle[1], eta_mle[2], eta_mle[3])
plot(tseq, msd_mle, type = "l")
points(tseq[1:1000], msd_fit(Xt, nlag = 1000, demean = TRUE),
       pch = 21, cex = .5, col = "red")


# ok now do this in a loop
Theta_mle <- sapply(unique(hbe$localid), function(ii) {
  message("Path ", ii)
  crouse_fit(ii)
})
Theta_mle <- t(Theta_mle)

# remove outliers
Theta_out <- Theta_mle[,-(1:2)]
Theta_out[,"tmax"] <- log10(Theta_out[,"tmax"])
iout <- rep(FALSE, nrow(Theta_out))
## iout <- Theta_out[,"tmax"] > 13
## iout <- iout | Theta_out[,"Sigma22"] > 5
## iout <- Theta_mle["Sigma22",] > 5

par(mfrow = c(3,3))
sapply(1:8, function(ii) {
  boxplot(Theta_out[!iout,ii], main = colnames(Theta_out)[ii])
})


# ok forget about outliers.
# Just simulate data from parameters as-is.

crouse_sim <- function(ii) {
  dT <- 1/60
  N <- Theta_mle[ii,"N"]
  eta <- Theta_mle[ii,c("alpha", "tmin", "tmax")]
  mu <- Theta_mle[ii,c("mu1", "mu2")]
  Sigma <- matrix(NA, 2, 2)
  Sigma[c(1,2,4)] <- Theta_mle[ii,c("Sigma11", "Sigma21", "Sigma22")]
  Sigma[1,2] <- Sigma[2,1]
  dX <- rSnorm(n = 2, acf = crouse_acf(eta, dT = dT, N = N)) %*% chol(Sigma)
  Xt <- rbind(0, apply(sweep(dX, 2, mu*dT, "+"), 2, cumsum))
  # round to instrument resolution
  res <- 0.00703125
  Xt <- round(Xt/res) * res
  data.frame(id = ii, x = Xt[,1], y = Xt[,2])
}

Xt <- crouse_sim(1)

hbe_sim <- do.call(rbind, lapply(unique(hbe$localid), crouse_sim))

# check whether data looks reasonable
hbe_sim %>% as_tibble %>%
  mutate(x = x + setNames(runif(n_distinct(id), -30, 30),
                          unique(id))[as.character(id)],
         y = y + setNames(runif(n_distinct(id), -30, 30),
                          unique(id))[as.character(id)]) %>%
  ggplot(aes(x = x, y = y, color = factor(id))) +
  geom_line()

# save
## hbe2 <- hbe
## hbe2$x <- hbe_sim$x
## write.csv(hbe_sim, file = "hbe_sim.csv", row.names = FALSE)
## write.csv(hbe2, file = "hbe2.csv", row.names = FALSE)
## hbe <- hbe_sim
## usethis::use_data(hbe)

hbe %>% as_tibble %>%
  group_by(localid) %>%
  summarize(xdiff = sort(unique(abs(diff(x))))[2],
            ydiff = sort(unique(abs(diff(y))))[2])
