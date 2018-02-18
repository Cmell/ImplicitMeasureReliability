setwd('~chrismellinger/GoogleDrive/ImplicitMeasureReliability/Sim3Results/')
load('GeneratedData0/DataIter110.RData')

findCor <- function (d) {
  d1 <- data.frame(snum=unique(d$snum))
  
  rts <- tapply(d$rt, INDEX=list(d$snum, d$half, d$pcat, d$tcat), mean, na.rm=T)
  browser()
  #d1 <- cast(d, snum ~ half + pcat + tcat, mean, value="rt")
  d1 <- data.frame(
    snum = dimnames(rts)[[1]]
  )
  d1 <- within(d1, {
    h1p1t1 = rts[snum,1,1,1]
    h1p1t2 = rts[snum,1,1,2]
    h1p2t1 = rts[snum,1,2,1]
    h1p2t2 = rts[snum,1,2,2]
    
    h2p1t1 = rts[snum,2,1,1]
    h2p1t2 = rts[snum,2,1,2]
    h2p2t1 = rts[snum,2,2,1]
    h2p2t2 = rts[snum,2,2,2]
    
    half1 <- (h1p1t1 - h1p1t2) - (h1p2t1 - h1p2t2)
    half2 <- (h2p1t1 - h2p1t2) - (h2p2t1 - h2p2t2)
  })
  rShT <- cor(d1$half1,d1$half2)
  return(rShT)
}

d$half <- 1
trialCount <- length(d$snum[d$snum==1])
sampledRows <- sample(1:trialCount, size=trialCount/2, replace=F)
for (s in unique(d$snum)) {
  d[d$snum==s, 'half'][sampledRows] <- 2
  #d[curRows, 'half'] <- rbinom(length(curRows), 1, .5)
}

findCor(d)
d$half <- ifelse(d$half==1, 2, 1)
findCor(d)

d$rnumx <- as.numeric(d$rnum) %% 2
d$rnumx[d$rnumx==0] <- 2
d1 <- cast(d, snum ~ rnumx + pcat + tcat, mean, value="rt")
d1$rnum1 <- (d1[,2] - d1[,3]) - (d1[,4] - d1[,5])
d1$rnum2 <- (d1[,6] - d1[,7]) - (d1[,8] - d1[,9])

# Build the return vector
r_pf <- cor(d1$rnum1,d1$rnum2)


# Find a -1 ====

flLst <- list.files(path='./', full.names = T, recursive = T)
flLst <- grep(paste0(resultDir, ".*csv"), flLst, value=T)

for (fl in flLst) {
  curRes <- read.csv(fl)
  if (curRes$r_sh < -.95) {
    print (fl)
    break
  }
}

# Graph some stuff ====
png(file = 'pHists.png', width=1600, height=1200)
op <- par(mfrow=c(5,6))

for (s in unique(d$snum)) {
  hist(d[d$snum==s, 'rt'], xlim = c(min(d$rt), max(d$rt)))
}
par(op)
dev.off()