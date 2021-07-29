#--- format the HBE data into a single csv file. -------------------------------

WT <- c("1p5", "2", "2p5", "3", "4", "5")

wt <- WT[1]
nm <- paste0("mucus_", wt, "_0.5um_")

datat <- as.matrix(read.csv(paste0(nm, "t.csv"))[,-1,drop=FALSE])
datax <- as.matrix(read.csv(paste0(nm, "x.csv"))[,-1,drop=FALSE])
datay <- as.matrix(read.csv(paste0(nm, "y.csv"))[,-1,drop=FALSE])

dataxy <- data.frame(id = rep(1:ncol(datay), each = nrow(datay)),
                     time = datat[,1], x = c(datax), y = c(datay))

id <- sample(ncol(datay), 1)
all(datax[,id] == dataxy[dataxy$id == id,"x"])
all(datay[,id] == dataxy[dataxy$id == id,"y"])

# ok bind all together
DataXY <- sapply(WT, function(wt) {
  nm <- paste0("mucus_", wt, "_0.5um_")
  datat <- as.matrix(read.csv(paste0(nm, "t.csv"))[,-1,drop=FALSE])
  datax <- as.matrix(read.csv(paste0(nm, "x.csv"))[,-1,drop=FALSE])
  datay <- as.matrix(read.csv(paste0(nm, "y.csv"))[,-1,drop=FALSE])
  data.frame(localid = rep(1:ncol(datay), each = nrow(datay)),
             wt = wt, time = datat[,1], x = c(datax), y = c(datay))
}, simplify = FALSE)

# global id
globalid <- sapply(WT, function(wt) {
  length(unique(DataXY[[wt]]$localid))
})
globalid <- cumsum(c(0, globalid[-length(globalid)]))
names(globalid) <- WT

# attach to DataXY and make into single data.frame
DataXY <- sapply(WT, function(wt) {
  cbind(globalid = globalid[wt] + DataXY[[wt]]$localid, DataXY[[wt]])
}, simplify = FALSE)
DataXY <- do.call(rbind, args = DataXY)

# save
# write.csv(DataXY, file = "HBE_data.csv", row.names = FALSE)

#--- compress to ship with package ---------------------------------------------

res <- 0.00000625 # instrument resolution
hbe <- read.csv("HBE_data.csv")
# this saves ~200Kb
hbe <- hbe[c("globalid", "wt", "x", "y")]
# this saves another ~200Kb
hbe$x <- round(hbe$x/res)*res
hbe$y <- round(hbe$y/res)*res

## save(hbe, file = "hbe.rda")
