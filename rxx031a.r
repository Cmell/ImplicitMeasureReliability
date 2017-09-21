# rxx03 series manipulates the variance of the prime distribution (a=0,b=1,c=2,d=3)

# Libraries ====

# add this to the search path for summit computing.
#.libPaths(c(.libPaths(), '/projects/chme2908/R_libs'))

library(CMUtils)

loadStuff(c(
  'car',
  #'psych',
  'reshape',
  #'gtheory',
  'parallel'
  #'pryr'
))

# Working Directory ====
setwd('~chrismellinger/GoogleDrive/ImplicitMeasureReliability/')
# setwd('/home/chme2908/ImplicitMeasureReliability/')

# Get arguments ====
library(optparse)

opts = 
  optionList = list(
    make_option(c("--nprim"), type="integer", help="number of primes"),
    make_option(c("--n.iter"), type="integer", help="number of iterations"),
    make_option(c("--ncores"), type="integer", help="number of cores to use"),
    make_option(c("--ntarg"), type="integer", help="number of targets"),
    make_option(c("--nsubj"), type="integer", help="number of subjects")
  )
optParser = OptionParser(option_list = optionList)
args = parse_args(optParser)

# Check the one required argument
if (is.null(args$n.iter)) {
  stop(
    paste("Must specify n.iter!")
  )
}
n.iter <- args$n.iter
# Default simulation parameters ====

# For the parrallelization
ncores <<- 2 * n.iter - 1 # This should be one less than is actually 
# available. Reserve one core for the main process.

# These are only needed if generating the data in files ahead of time.
# At present, this method is not used.
# scratchDirLo <- './scratch/LowVarData'
# scratchDirHi <- './scratch/HighVarData'

nsubj <<- 32
nprim <<- 4
npcat <<- 2
ntarg <<- 4
ntcat <<- 2
nreps <<- 2 # should be even number

svar <<- 1
pvarLo <<- 1
pvarHi <<- pvarLo * 2
tvar <<- 1
evar <<- 1

# Substitute Provided Arguments for Default Values ====

# Overwrite the defaults when they are provided.
for (var in names(args)) {
  assign(var, args[var][[1]])
}

# Print Important Parameters for the Logs ====
varLst <- c(
  "ncores",
  "n.iter",
  "nsubj",
  "nprim",
  "npcat",
  "ntarg",
  "ntcat",
  "nreps",
  "svar",
  "pvarLo",
  "pvarHi",
  "tvar",
  "evar"
)
for (var in varLst) {
  print(paste0(var, ": ", get(var)))
}

# Random Seed ====

# Random number considerations. This will make the result reproducible 
# and also ensure that each iteration is reasonably independent.
RNGkind("L'Ecuyer-CMRG")
set.seed(593065038)

# seed. Deprecated. This is now accomplished when creating the workers.
# set.seed(12345)   

# Date String ====

dateStr <- format(Sys.time(), format='%Y-%m-%d_%H.%M.%S')
startTime <<- proc.time()["elapsed"]
# profileFl <- paste0(dateStr, '_profile.txt')

# Build functions for the job ====

genData = function(
                  nsubj,
                  nprim,
                  npcat,
                  ntarg,
                  ntcat,
                  nreps,# should be even number

                  svar,
                  pvar,
                  tvar,
                  evar
                  ) {
  # subject differences ====
  
  snum <<- as.factor(1:nsubj)
  pnum <<- as.factor(1:(npcat*nprim))
  tnum <<- as.factor(1:(ntcat*ntarg))
  rnum <<- as.factor(1:nreps)
  
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
  d <- expand.grid(rnum=rnum,tnum=tnum,pnum=pnum,snum=snum)
  d <- merge(d,target,by.x="tnum",by.y="tnum")
  d <- merge(d,prime,by.x="pnum",by.y="pnum")
  d <- merge(d,subj,by.x="snum",by.y="snum")
  
  
  d$error <- rnorm(nrow(d),0,evar)
  # d$rt <- 600 + 1*(d$pcat*d$tcat) + 1*(d$pcat*d$tcat*d$prej) + 
  #   1*(d$pcat*d$tcat*d$prej*d$tvaln) + 1*(d$pcat*d$tcat*d$prej*d$pprot) + 
  #   1*(d$pcat*d$tcat*d$prej*d$pprot*d$tvaln) + 5*d$error
  
  d$rt <- 6.4 + 0.07*(d$prej) + 0.00*(d$pcat) + 0.02*(d$tcat) + 0.02*(d$pprot) + 
    0.02*(d$tvaln) + 0.00*(d$prej*d$pcat) + 0.00*(d$prej*d$tcat) + 
    0.04*(d$pcat*d$tcat) + 0.04*(d$prej*d$pprot) + 0.02*(d$pprot*d$tcat) + 
    0.00*(d$prej*d$tvaln) + 0.00*(d$pcat*d$tvaln) + 0.00*(d$pprot*d$tvaln) + 
    0.08*(d$prej*d$pcat*d$tcat) + 0.00*(d$prej*d$pcat*d$tcat*d$tvaln) + 
    0.06*(d$prej*d$pcat*d$tcat*d$pprot) +
    0.05*(d$pcat*d$tcat*d$prej*d$pprot*d$tvaln) + 0.25*d$error
  
  
  # variance decomposition ====
  d$snum <- as.factor(d$snum)
  d$pcat <- as.factor(d$pcat)
  d$pnum <- as.factor(d$pnum)
  d$tcat <- as.factor(d$tcat)
  d$tnum <- as.factor(d$tnum)
  
  # Why are pcat & tcat factors? Looks like they are randomly generated.
  
  # contrasts(d$snum) <- contr.poly
  # contrasts(d$pcat) <- contr.poly
  # contrasts(d$pnum) <- contr.poly
  # contrasts(d$tcat) <- contr.poly
  # contrasts(d$tnum) <- contr.poly
  
  return(d)
}

modelFn <- function (d, i=-1) 
{
  m1 <- lm(rt ~ snum * pcat * tcat + snum * pnum * tcat + 
             snum * pcat * tnum + 
             snum * pnum * tnum, data=d)
  print(paste0('Iteration ', i, ' lm() call finished.'))
  
  my.anova1 <- Anova(m1, type="III", singular.ok = TRUE)
  print(paste0('Iteration ', i, ' Anova() call finished.'))
  numPar <- nrow(my.anova1)
  est <- my.anova1[1:numPar,'Sum Sq'] / my.anova1[1:numPar,'Df']
  names(est) <- rownames(my.anova1)
  
  # split-half reliability

  d$half <- 1
  d$half[sample(1:nrow(d), size=nrow(d)/2, replace=F)] <- 2
  
  d1 <- cast(d, snum ~ half + pcat + tcat, mean, value="rt")
  d1$half1 <- (d1[,2] - d1[,3]) - (d1[,4] - d1[,5])
  d1$half2 <- (d1[,6] - d1[,7]) - (d1[,8] - d1[,9])
  est['r_sh'] <- cor(d1$half1,d1$half2)

  # parallel forms reliability, must recode replications to only 
  # 2 values (odd and even)
  d$rnumx <- as.numeric(d$rnum) %% 2
  d$rnumx[d$rnumx==0] <- 2
  d1 <- cast(d, snum ~ rnumx + pcat + tcat, mean, value="rt")
  d1$rnum1 <- (d1[,2] - d1[,3]) - (d1[,4] - d1[,5])
  d1$rnum2 <- (d1[,6] - d1[,7]) - (d1[,8] - d1[,9])
  
  # Build the return vector
  est['r_pf'] <- cor(d1$rnum1,d1$rnum2)
  
  rm(modForm, m1, my.anova1, numPar, d, d1)
  return(est)
}

finishCalcs <- function (est) {
  est$var.snpctc <- 
    (est$ms.snpctc - est$ms.snpctn - est$ms.snpntc + est$ms.snpntn) / 
    (nprim*ntarg*nreps)
  est$var.snpctn <- (est$ms.snpctn - est$ms.snpntn) / (nprim*nreps)
  est$var.snpntc <- (est$ms.snpntc - est$ms.snpntn) / (ntarg*nreps)
  est$var.snpntn <- (est$ms.snpntn - est$ms.resid) / (nreps)
  est$var.resid <- est$ms.resid
  
  # recode negative variances to zero (this will bias the estimates but is 
  # necessary to avoid negative reliabilities)
  est$var.snpctc[est$var.snpctc<0] <- 0
  est$var.snpctn[est$var.snpctn<0] <- 0
  est$var.snpntc[est$var.snpntc<0] <- 0
  est$var.snpntn[est$var.snpntn<0] <- 0
  
  est$rxxmse <- (est$ms.snpctc-est$ms.resid) / est$ms.snpctc
  
  est$rxxvar <- 
    (est$var.snpctc + (est$var.snpctn / ntarg) + (est$var.snpntc/ nprim) + 
       (est$var.snpntn / (nprim*ntarg))
    ) / 
    (est$var.snpctc + (est$var.snpctn / ntarg) + (est$var.snpntc / nprim) + 
       (est$var.snpntn / (nprim*ntarg)) + est$var.resid / (nprim*ntarg*nreps)
    )
  
  est$rxxvar.prand <- 
    est$var.snpctc / 
    (est$var.snpctc + (est$var.snpctn / ntarg) + (est$var.snpntc / nprim) + 
       (est$var.snpntn / (nprim*ntarg)) + est$var.resid / (nprim*ntarg*nreps)
    )
  return(est)
}

renameEstCols <- function (est) {
  colnames(est) <- c(
    "ms.int","ms.sn","ms.pc","ms.tc","ms.pn","ms.tn","ms.snpc","ms.sntc",
    "ms.pctc","ms.snpn","ms.tcpn","ms.sntn","ms.pctn","ms.pntn","ms.snpctc",
    "ms.snpntc","ms.snpctn","ms.snpntn","ms.resid","r_sh","r_pf", "var"
    )
  return(est)
}

iterFn <- function (i, curPvar) {
  print(paste('Beginning iteration', i, sep=" "))
  
  initTm <- proc.time()[3]
  d <- genData(
    nsubj = nsubj,
    nprim = nprim,
    npcat = npcat,
    ntarg = ntarg,
    ntcat = ntcat,
    nreps = nreps, # should be even number
    
    svar = svar,
    pvar = curPvar,
    tvar = tvar,
    evar = evar
  )
  tm <- proc.time()[3] - initTm
  print(paste0('Iteration ', i, ' data gen time: ', tm))
  
  initTm <- proc.time()[3]
  curEst <- modelFn(d, i=i)
  tm <- proc.time()[3] - initTm
  #print(paste0('Iteration ', i, ' model time: ', tm))
  
  rm(d, initTm, tm); gc();
  return(curEst)
}

 
# # Parallelize! the modeling part... ====

mc.reset.stream()

# This one does the low variance simulation.
# Rprof(profileFl)
system.time({
  estLst <- mclapply(1:(2 * n.iter), 
                     function (i) {
                       # If i is even, do the low ones, if it is odd, the high
                       if (i %% 2 == 0) {
                         return(
                           c(iterFn(i, curPvar=pvarLo), var=pvarLo)
                         )
                       } else {
                         return(
                           c(iterFn(i, curPvar=pvarHi), var=pvarHi)
                         )
                       }
                     },
                     mc.cores = ncores
  )
})
# Rprof(NULL)
# summaryRprof(profileFl)

# Make the est dataframe. ====

est <- data.frame(Reduce(rbind, estLst))

est <- renameEstCols(est)

# Perform calculations. ====

est <- finishCalcs(est)

# Save all the things. ====

print(paste('Run at', dateStr))

save(est,
     file=paste0('est_', dateStr,'.RData')
     )

# Testing ====

# d <- genData(
#   nsubj = nsubj,
#   nprim = nprim,
#   npcat = npcat,
#   ntarg = ntarg,
#   ntcat = ntcat,
#   nreps = nreps, # should be even number
#   
#   svar = svar,
#   pvar = 1,
#   tvar = tvar,
#   evar = evar
# )

# Using clusters - deprecated ====

# if (ncores > n.iter) {
#   numCores <- n.iter
# } else {
#   numCores <- ncores
# }

# cl <- makeCluster(numCores, outfile=paste0(dateStr, '.out'))
# clusterExport(cl, ls())
# clusterEvalQ(cl, {
#   #.libPaths(c(.libPaths(), '/projects/chme2908/R_libs'))
#   library(car)
#   library(reshape)
#   library(psych)
# })