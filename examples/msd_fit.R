# Compute the MSD of Xt
nlag <- 600
msd <- msd_fit(Xt, nlag = nlag)

plot(dt*(1:nlag), msd, xlab = "Time", ylab = "MSD")

