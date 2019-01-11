require(subdiff)
source(system.file("proj", "msd-functions.R", package = "subdiff"))
source(system.file("proj", "drift_removal_func.R", package = "subdiff"))

# specify the absolute path.
filePath <- "D:/Dropbox/gpsi/JOR/5MD PEO and sucrose dec 12 2018"
samples <- dir(file.path(filePath, "rds files"))

# for example we are interested in PEO5MD 0p6
dataID <- 9
paths <- readRDS(file.path(filePath, "rds files",samples[dataID]))

# movie information and filter
movie_info <- movie.info(paths)
movie_seq <- sort(unique(movie_info$info))
movie_count <- rep(NA, length(movie_seq))
for(ii in 1:length(movie_seq)) {
  movie_count[ii] <- sum(movie_info$info == movie_seq[ii])
}
movie_count <- rbind(movie_seq, movie_count)
rownames(movie_count) <- c("Movie ID", "# of good paths")
# show the movie IDs and number of good paths for each ID
movie_count


# for one movie -----------------------------------------------------------

# extract the movie-wise drift for movie 8
movieID <- 8
movie_drift <- movie.drift(paths, movieID, movie_info)
drift <- apply(movie_drift$drift, 2, cumsum)
Xt <- apply(movie_drift$incr_x, 2, cumsum)
Yt <- apply(movie_drift$incr_y, 2, cumsum)

# plot the trajectories
xrng <- range(Xt)
yrng <- range(Yt)
ll <- ncol(Xt)
clrs <- rainbow(ll)
plot(x = 0, type = "n", xlim = xrng, ylim = yrng, xlab = "x-axis", ylab = "y-axis", 
     main = "movie 1, 1p52")
for(ii in 1:ll) {
  lines(x = Xt[,ii], y = Yt[,ii], col = clrs[ii])
}

# plot the movie drift
plot(x = drift[,1], y = drift[,2], type = "l", xlab = "", ylab = "")

# plot the de-trended trajectories
xrng <- range(Xt - drift[,1])
yrng <- range(Yt - drift[,2])
plot(x = 0, type = "n", xlim = xrng, ylim = yrng, xlab = "x-axis", ylab = "y-axis", 
     main = "movie 1, 1p52")
for(ii in 1:ll) {
  lines(x = Xt[,ii]-drift[,1], y = Yt[,ii]-drift[,2], col = clrs[ii])
}


# estimation:

# original
est0 <- movie.est(paths, movieID, movie_info, movie_drift, NL.method = FALSE)
# drift-removed
est1 <- movie.est(paths, movieID, movie_info, movie_drift)

# plot the estimation results
xrng <- range(est1[,1], est0[,1])
yrng <- range(est1[,2], est0[,2])
lbl <- strsplit(samples[dataID], '5MD')[[1]][2]
lbl <- strsplit(lbl, 'pct')[[1]][1]
plot(x = 1, type = "n", xlim = xrng, ylim = yrng, xlab = "alpha", ylab = "D", 
     main = paste0("movie ", movieID, ", ", lbl),log = "y")
points(x = est0[,1], y = est0[,2], col = clrs, pch = 19)
points(x = est1[,1], y = est1[,2], col = clrs, pch = 15)
segments(x0 = est1[,1], x1 = est0[,1],
         y0 = est1[,2], y1 = est0[,2], lty = 2,
         col = rgb(t(col2rgb(clrs)), alpha = 100, maxColorValue = 255))
legend("topleft", legend = c("original", "de-trend"), pch = c(19,15), 
       bty = "n", cex = 0.6)

# it is possible that some path (the yellow one in this movie) does not quite follow the 
# movie-wise drift.
movie_R2 <- movie.R2(movie_drift)
# R^2 for this movie
median(movie_R2)

# for each path, label its R^2 on the plot
text(xy.coords(x = est0[,1], y = est0[,2]), labels = round(movie_R2,3))

# So far what we obtain is the movie-wise estimation, if we want the data-wise estimation, 
# we first compute this for all movies, and use the median/mean of estimation for this data
# when number of paths in one data is small, median is more robust than mean.

# since 1p52 only has one movie (movie 1) with more than 5 paths, the estimation of 1p52 is
apply(est1, 2, median)
# and the original method:
apply(est0, 2, median)


# for whole data -----------------------------------------------------------

# find all movies with at least 5 good paths
movie_seq <- (movie_count[1,])[movie_count[2,]>4]

# compute estimation results for all movies
est_old <- est_new <- numeric()
for(movieID in movie_seq) {
  # progress bar
  print(paste0("Movie ID: ", movieID))
  movie_drift <- movie.drift(paths, movieID, movie_info)
  
  est0 <- movie.est(paths, movieID, movie_info, movie_drift, NL.method = FALSE)
  est1 <- movie.est(paths, movieID, movie_info, movie_drift)

  # decide using est0 or est1 by checking the value of R^2
  # is 0.2 a proper value?
  movie_R2 <- movie.R2(movie_drift)
  if(median(movie_R2) < 0.2) {
    est_new <- rbind(est_new, est0)
  } else {
    est_new <- rbind(est_new, est1)
  }
  est_old <- rbind(est_old, est0)
}

xrng <- range(est_new[,1], est_old[,1])
yrng <- range(est_new[,2], est_old[,2])
lbl <- strsplit(samples[dataID], '5MD')[[1]][2]
lbl <- strsplit(lbl, 'pct')[[1]][1]
clrs <- rainbow(nrow(est_new))
plot(x = 1, type = "n", xlim = xrng, ylim = yrng, xlab = "alpha", ylab = "D", 
     main = paste0(lbl),log = "y")
points(x = est_old[,1], y = est_old[,2], col = clrs, pch = 19)
points(x = est_new[,1], y = est_new[,2], col = clrs, pch = 15)
segments(x0 = est_new[,1], x1 = est_old[,1],
         y0 = est_new[,2], y1 = est_old[,2],
         lty = 2,
         col = rgb(t(col2rgb(clrs)), alpha = 100, maxColorValue = 255))
legend("topleft", legend = c("original", "de-trend"), pch = c(19,15), 
       bty = "n", cex = 0.6)

# data-wise estimation
apply(est_new, 2, mean)

