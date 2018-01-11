require(subdiff)
require(doParallel) # parallelize

data_path <- "C:/Users/Yun/Dropbox/gpsi/Inference_Results_for_Viscous_Data/Data"
data_name <- "Water60"
# load statistics
print(load(file.path(data_path, "model-stats",
                     paste0("fma-stats_", data_name, ".RData"))))

#-------------------------------------------------------------------------------
# calculate the point and interval estimate of alpha_trans and D_trans
npaths <- length(data_stats)
N <- rep(NA, npaths)
movieID <- N
for(ii in 1:npaths) {
  ds <- data_stats[[ii]]
  if(!anyNA(ds)) {
    N[ii] <- ds$size
    movieID[ii] <- ds$movieID
  }
}
# remove the NA terms
filt <- !is.na(N) & N > 1000
N <- N[filt]
movieID <- movieID[filt]
npaths <- length(N)
nmovies <- length(unique(movieID))
# parameter for simulation
alpha <- 1
rho <- 0.3
dT <- 1/60

fma_sim <- function(alpha, rho, dT, N) {
  dX <- rSnorm(n = 2, acf = fbm_acf(alpha, dT, N))
  dY <- dX * (1-rho) + rbind(0, dX[-N, ]) * rho
  dY
}

fma_data <- sapply(1:npaths, function(ii) {
  fma_sim(alpha, rho, dT, N[ii])
})

# inference
nclust <- detectCores(logical = FALSE)
cl <- makeCluster(nclust)
registerDoParallel(cl)
clusterSetRNGStream(cl = cl)
cl_pkg <- c("subdiff", "nortest")

system.time({
  data_stats <- foreach(ii = 1:npaths, .packages = cl_pkg) %dopar% {
    dX <- fma_data[[ii]]
    fit <- try(fma_fit(dX = dX, dT = dT), silent = TRUE) # fit the model
    list(coef = fit$coef, vcov = fit$vcov)
  }
})

stopCluster(cl)

# cochran Q test accross the movies
#-------------------------------------------------------------------------------
# calculate the point and interval estimate of alpha_trans and D_trans
ncoef <- 7
coef <- matrix(NA, npaths, ncoef)
vcov <- array(NA, c(ncoef, ncoef, npaths))

for(ii in 1:npaths) {
  ds <- data_stats[[ii]]
  if(!anyNA(ds)) {
    coef[ii,] <- ds$coef
    vcov[,,ii] <- ds$vcov
  }
}

# cochran Q test ----------------------------------------------------------
cochran_stat <- function(ind) {
  filt <- movieID == ind
  coef_avg <- coef[filt, ]
  vcov_avg <- vcov[, , filt]
  stat <- cochranMQ(coef_avg, vcov_avg)
  list(coef = stat$est, vcov = stat$ve)
}

movie_coef <- matrix(NA, nmovies, ncoef)
movie_vcov <- array(NA, c(ncoef, ncoef, nmovies))

for(ii in 1:nmovies) {
  print(ii)
  ds <- cochran_stat(ii)
  movie_coef[ii,] <- ds$coef
  movie_vcov[,,ii] <- ds$vcov
}

Q_stat <- cochranMQ(movie_coef, movie_vcov)
Q_stat$pval

# movie-level estimation of alpha-D plot
plot(x = movie_coef[, 1], y = movie_coef[, 5], xlab = expression(logit(alpha/2)), 
     ylab = expression(log(D)), col = "blue")
# original scale
plot(x = itrans_alpha(movie_coef[, 1]), y = exp(movie_coef[, 5]), xlab = expression(logit(alpha/2)), 
     ylab = expression(log(D)), col = "blue")
# path-wise estimation of alpha-D plot
plot(x = coef[, 1], y = coef[, 5], xlab = expression(logit(alpha/2)), 
     ylab = expression(log(D)), col = "blue")
