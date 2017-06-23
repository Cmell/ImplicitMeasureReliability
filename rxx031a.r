# rxx03 series manipulates the variance of the prime distribution (a=0,b=1,c=2,d=3)
# 32 subjects, 4 primes, 4 targets, 2 reps, pvar=1

# libraries ====
if(!require(devtools)) {
  install.packages('devtools')
  library(devtools)
}
install_github('Cmell/CMUtils')
library(CMUtils)

loadStuff(c(
  'car',
  'psych',
  'reshape',
  'gtheory',
  'parallel'
))

# seed
# set.seed(12345)   

iterFn <- function (i) 
{
  
  # subject differences ====
  # This is fast.
  prej <- rnorm(snum,0,svar)
  basert <- rnorm(snum,0,1)
  subj <- cbind(snum,prej,basert)
  
  
  # prime differences ====
  # This is fast.
  pprot <- rnorm((npcat*nprim),0,pvar)
  pcat <- rep(rnorm(npcat,0,1),nprim)
  prime <- cbind(pnum,pcat,pprot)
  
  
  # target differences ====
  # this is fast
  tvaln <- rnorm((ntcat*ntarg),0,tvar)
  tcat <- rep(rnorm(ntcat,0,1),ntarg)
  target <- cbind(tnum,tcat,tvaln)
  
  
  # build integrated data file ====
  # This is nontrivially slower, but still pretty fast.
  d <- expand.grid(rnum=rnum,tnum=tnum,pnum=pnum,snum=snum)
  d <- merge(d,target,by.x="tnum",by.y="tnum")
  d <- merge(d,prime,by.x="pnum",by.y="pnum")
  d <- merge(d,subj,by.x="snum",by.y="snum")
  
  
  # This is fast.
  d$error <- rnorm(nrow(d),0,evar)
  d$rt <- 600 + 1*(d$pcat*d$tcat) + 1*(d$pcat*d$tcat*d$prej) + 
    1*(d$pcat*d$tcat*d$prej*d$tvaln) + 1*(d$pcat*d$tcat*d$prej*d$pprot) + 
    1*(d$pcat*d$tcat*d$prej*d$pprot*d$tvaln) + 5*d$error
  

  # variance decomposition ====
  d$snum <- as.factor(d$snum)
  d$pcat <- as.factor(d$pcat)
  d$pnum <- as.factor(d$pnum)
  d$tcat <- as.factor(d$tcat)
  d$tnum <- as.factor(d$tnum)
  
  contrasts(d$snum) <- contr.poly
  contrasts(d$pcat) <- contr.poly
  contrasts(d$pnum) <- contr.poly
  contrasts(d$tcat) <- contr.poly
  contrasts(d$tnum) <- contr.poly

  # This step is worth about 46 seconds.
  m1 <- lm(rt ~ snum * pcat * tcat + snum * pnum * tcat + snum * pcat * tnum + 
             snum * pnum * tnum, data=d)
  
  # This step is worth about 63 seconds
  my.anova1 <- Anova(m1, type="III", singular.ok = TRUE)
  numPar <- nrow(my.anova1)
  est <- my.anova1[1:numPar,'Sum Sq']/my.anova1[1:numPar,'Df']
  names(est) <- rownames(my.anova1)
  
  
  # split-half reliability

  d$half <- 1
  d$half[sample(1:nrow(d), size=nrow(d)/2, replace=F)] <- 2
  
  d1 <- cast(d, snum ~ half + pcat + tcat, mean, value="rt")
  d1$half1 <- (d1[,2] - d1[,3]) - (d1[,4] - d1[,5])
  d1$half2 <- (d1[,6] - d1[,7]) - (d1[,8] - d1[,9])
  est['r_sh'] <- cor(d1$half1,d1$half2)

  # parallel forms reliability, must recode replications to only 2 values (odd and even)
  d$rnumx <- as.numeric(d$rnum) %% 2
  d$rnumx[d$rnumx==0] <- 2
  d1 <- cast(d, snum ~ rnumx + pcat + tcat, mean, value="rt")
  d1$rnum1 <- (d1[,2] - d1[,3]) - (d1[,4] - d1[,5])
  d1$rnum2 <- (d1[,6] - d1[,7]) - (d1[,8] - d1[,9])
  
  # Build the return vector
  est['r_pf'] <- cor(d1$rnum1,d1$rnum2)
  
  return(est)
}

# set numbers ====
n.iter <<- 6
nsubj <<- 32
nprim <<- 4
npcat <<- 2
ntarg <<- 4
ntcat <<- 2
nreps <<- 2 # should be even number

svar <<- 1
pvar <<- 1
tvar <<- 1
evar <<- 1

snum <<- as.factor(1:nsubj)
pnum <<- as.factor(1:(npcat*nprim))
tnum <<- as.factor(1:(ntcat*ntarg))
rnum <<- as.factor(1:nreps)


# Simple function call ====
#system.time(estLst <- lapply(1:n.iter, iterFn))
# 211 seconds elapsed.

# Parallelize! ====
numCoresAvail <- detectCores()
if (numCoresAvail > n.iter) {
  numCores <- n.iter
} else {
  numCores <- numCoresAvail
}

cl <- makeCluster(numCores)
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(car)
  library(reshape)
  library(psych)
})
system.time({
  estLst <- parLapply(cl, 1:n.iter, iterFn)
})
stopCluster(cl)
# Nice! 109 seconds elapsed. Basically half the time, which is what should happen.

# Later math ====

# Make the est dataframe.
est <- data.frame(Reduce(rbind, estLst))

colnames(est) <- c("ms.int","ms.sn","ms.pc","ms.tc","ms.pn","ms.tn","ms.snpc","ms.sntc","ms.pctc","ms.snpn","ms.tcpn","ms.sntn","ms.pctn","ms.pntn","ms.snpctc","ms.snpntc","ms.snpctn","ms.snpntn","ms.resid","r_sh","r_pf")

est$var.snpctc <- (est$ms.snpctc - est$ms.snpctn - est$ms.snpntc + est$ms.snpntn)/(nprim*ntarg*nreps)
est$var.snpctn <- (est$ms.snpctn - est$ms.snpntn)/(nprim*nreps)
est$var.snpntc <- (est$ms.snpntc - est$ms.snpntn)/(ntarg*nreps)
est$var.snpntn <- (est$ms.snpntn - est$ms.resid)/(nreps)
est$var.resid <- est$ms.resid

# recode negative variances to zero (this will bias the estimates but is necessary to avoid negative reliabilities)
est$var.snpctc[est$var.snpctc<0] <- 0
est$var.snpctn[est$var.snpctn<0] <- 0
est$var.snpntc[est$var.snpntc<0] <- 0
est$var.snpntn[est$var.snpntn<0] <- 0

est$rxxmse <- (est$ms.snpctc-est$ms.resid)/est$ms.snpctc
est$rxxvar <- (est$var.snpctc+(est$var.snpctn/ntarg)+(est$var.snpntc/nprim)+(est$var.snpntn/(nprim*ntarg)))/(est$var.snpctc+(est$var.snpctn/ntarg)+(est$var.snpntc/nprim)+(est$var.snpntn/(nprim*ntarg))+est$var.resid/(nprim*ntarg*nreps))
est$rxxvar.prand <- est$var.snpctc/(est$var.snpctc+(est$var.snpctn/ntarg)+(est$var.snpntc/nprim)+(est$var.snpntn/(nprim*ntarg))+est$var.resid/(nprim*ntarg*nreps))
describe(est)

save(est, file='est.RData')
#est1 <- rbind(est1, est)
#write.csv(est1,"rxx031_est.csv")
