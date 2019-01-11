# extract the filter and movie information from paths
# filter: removing SNR > 4 and length < 1145 (the max length)
# Is it necessary to remove all paths shorter than 1145?
movie.info <- function(paths) {
  pathID <- unique(paths$pathID)
  filter_p <- numeric()
  N <- 1145
  for (jj in 1:length(pathID)) {
    index <- paths$pathID == pathID[jj]
    SNR <- paths$snr[index]
    # x_incr <- diff(paths$X[index])
    # y_incr <- diff(paths$Y[index])
    if (all(SNR > 4) && length(SNR) == N) {
      filter_p <- c(filter_p, pathID[jj])
    }
  }
  npaths <- length(filter_p)
  # get the movie information
  movie_seq <- rep(NA, npaths)
  for (jj in 1:npaths) {
    index <- paths$pathID == filter_p[jj]
    movie_seq[jj] <- unique(paths$movieID[index])
  }
  list(info = movie_seq, filter = filter_p)
}

# functions for extracting the movie-wise mean for specific movie
# movie_ID is in movie_seq
# incr_drift is only the increment of movie-drift
movie.drift <- function(paths, movie_ID, movie_info) {
  N <- 1145
  
  if(missing(movie_info)) movie_info <- movie.info(paths)
  
  movie_seq <- movie_info$info
  filter_p <- movie_info$filter
  index <- filter_p[movie_seq == movie_ID]
  
  # extract the x-axis and y-axis increments
  # only for movies with more than 5 paths
  ll <- length(index)
  if(ll > 5) {
    dX <- dY <- matrix(NA, N-1, ll)
    incr_drift <- matrix(NA, N-1, 2)
    for(kk in 1:ll) {
      incr <- diff(get_path(index[kk], paths = paths))
      dX[,kk] <- incr[,1]     
      dY[,kk] <- incr[,2]
    }
    incr_drift[,1] <- apply(dX, 1, mean)
    incr_drift[,2] <- apply(dY, 1, mean)
    
    list(incr_x = dX, incr_y = dY, drift = incr_drift)
  } else {
    stop("This movie does not have sufficient good paths ")
  }
}

# estimating subdiff-parameters (alpha, D) for specific movie
# using fma model as default
# if NL.method, using movie-wise drift removal, otherwise linear drift removal
movie.est <- function(paths, movie_ID, movie_info, movie_drift, NL.method = TRUE) {
  N <- 1145
  dT <- 1/38.17
  
  if(missing(movie_info)) movie_info <- movie.info(paths)
  if(missing(movie_drift)) movie_drift <- movie.drift(paths, movie_ID, movie_info)
  dX <- movie_drift$incr_x
  dY <- movie_drift$incr_y
  ll <- ncol(dX)
  
  # removing the movie-wise drift
  if(NL.method) {
    dX <- dX - movie_drift$drift[,1]
    dY <- dY - movie_drift$drift[,2]
  }
  
  est <- matrix(NA, ll, 2)
  colnames(est) <- c("alpha", "D")
  for(ii in 1:ll) {
    dZ <- cbind(dX[,ii], dY[,ii])
    ans <- fma_fit(dZ, dT, var_calc = FALSE)
    est[ii, 1] <- itrans_alpha(ans["gamma"])
    est[ii, 2] <- exp(ans["lambda1"])
  }
  est
}

# computation of R^2:
# particle = drift + subdiffusion, 
# var(particle) = E[var(particle | drift)] + var(E[particle | drift])
# E[var(particle | drift)] = E[var(subdiffusion)]
# R^2 = 1 - E[var(subdiffusion)] / var(particle)
movie.R2 <- function(movie_drift) {
  dX <- movie_drift$incr_x
  dY <- movie_drift$incr_y
  incr_drift <- movie_drift$drift
  ll <- ncol(dX)
  
  v1 <- var(c(dX, dY))
  e1 <- rep(NA, ll)
  for(ii in 1:ll) {
    e1[ii] <- var(c(dX[,ii] - incr_drift[,1], dY[,ii] - incr_drift[,2]))
  }
  1 - e1 / v1
}
