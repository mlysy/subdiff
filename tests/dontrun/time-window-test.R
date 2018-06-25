# plot the estimation
require(subdiff)

data_path <- "D:/Dropbox/gpsi/cmerr-outline"
source(paste0(data_path, "/R/plot_func.R"))

show_name <- expression(tau[0] == 1/500, tau[0] == 1/270,
                        tau[0] == 1/150, tau[0] == 1/80,
                        tau[0] == 1/42, tau[0] == 1/23, tau[0] == 1/13)


# gle, plot ---------------------------------------------------------------

# compare the share of GLE
dT <- 1/60
N <- 1800
alpha <- 0.5
tau <- exp(seq(log(1e-4), log(8e-2), len = 7))
K <- 300
nsize <- length(tau)

tSeq <- function(lim, N = 1800, log = TRUE) {
  if(log) lim <- log(lim)
  ts <- seq(lim[1], lim[2], length.out = N)
  if(log) ts <- exp(ts)
  ts
}

# MSD plot
msd_mat1 <- matrix(NA, N, nsize)
tfseq <- range(rouse_sub(alpha, tau[1], K)[1]/10, dT, 
               rouse_sub(alpha, tau[1], K)[2]*10, N*dT)
tfseq <- tSeq(tfseq, N)
for(jj in 1:nsize) {
  msd_mat1[,jj] <- rouse_msd(tfseq, alpha, tau[jj], K)
}

clrs <- rainbow(nsize)
par(mar = c(1,1,3,0)+.1)
multi.plot(x = tfseq, y = msd_mat1, clrs = clrs, 
           xlab = "Time Sequence", ylab = "MSD", xaxt = "n")
for(jj in 1:nsize) {
  abline(v = rouse_sub(alpha, tau[jj], K)[1:2], col = clrs[jj], lty = 2)
}
abline(v = c(1/60, 30), col = "black", lty = 3)


# on normal scale
tseq <- 1:N * dT
tmin <- rep(NA, nsize)
alpha1 <- rep(NA, nsize)
D1 <- rep(NA, nsize)
msd_mat2 <- matrix(NA, N, 1)
for(jj in 1:1) {
  msd_mat2[,jj] <- rouse_msd(tseq, alpha, tau[jj], K)
}
theta <- time_window(msd_mat2, tseq, error = 0.03, tmax = TRUE, log = TRUE)
theta2 <- time_window(msd_mat2, tseq, error = 0.03, tmax = TRUE)

par(mfrow = c(1,1), mar = c(3,3,2,0)+.1, mgp = c(1.5,0.5,0), oma=c(0,0,1,0))


multi.plot(x = tfseq, y = as.matrix(msd_mat1[,1]), clrs = clrs, 
           xlab = "Time Sequence", ylab = "MSD", 
           main = "MSD of GLE simulation")
par(new = TRUE)
plot(x = tfseq, y = theta2[4]*tfseq^theta2[3], col = "blue", 
     log = "xy", lty = 2, ylim = range(msd_mat1[,1]), type = "l", 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
abline(v = theta2[1:2], col = "blue", lty = 3)
par(new = TRUE)
plot(x = tfseq, y = theta[4]*tfseq^theta[3], col = "purple", 
     log = "xy", lty = 2, ylim = range(msd_mat1[,1]), type = "l", 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
abline(v = theta[1:2], col = "purple", lty = 3)
abline(v = c(1/60, 30), col = "black", lty = 2)
legend("bottomright", col = c("blue", "purple"), 
       legend = c("normal scale", "log scale"), lty = 2, bty = "n")
