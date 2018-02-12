# rxx03 series manipulates the variance of the prime distribution (a=0,b=1,c=2,d=3)

# Libraries ====

# add this to the search path for summit computing.
#.libPaths(c(.libPaths(), '/projects/chme2908/R_libs'))

#library(CMUtils)
pkgLst <- c(
  'car',
  'optparse',
  #'psych',
  'reshape',
  'pbdMPI'
  #'filelock'
  #'gtheory',
  #'parallel'
  #'pryr'
)
for (p in pkgLst) {
  library(p, character.only = T, quietly = T)
}

# Date String ====

suppressWarnings({
  dateStr <- format(Sys.time(), format='%Y-%m-%d_%H.%M.%S')
  comm.print(paste('Run start time:', format(Sys.time(), format='%Y-%m-%d_%H:%M:%S')))
  startTime <<- proc.time()["elapsed"]
})

# Working Directory ====
cmDir <- '~chrismellinger/GoogleDrive/ImplicitMeasureReliability/'
corcDir <- '/scratch/summit/chme2908/ImplicitMeasureReliability/'
if (dir.exists(cmDir)) {
  setwd(cmDir)
} else if (dir.exists(corcDir)) {
  setwd(corcDir)
}

# Get arguments ====

opts = 
  optionList = list(
    make_option(c("--nprim"), type="integer", help="number of primes"),
    make_option(c("--n.iter"), type="integer", 
                help="
                Number of iterations (the iteration number to end on). Note that
                this is multiplied by the number of variances in the variance
                list for the final number of iterations.
                "),
    make_option(c("--iter.start"), type="integer", 
                help="
                Iteration number to begin at. This is with reference to the
                total number of iterations n.iter * length(varianceLst).
                "),
    make_option(c("--ncores"), type="integer", help="number of cores to use"),
    make_option(c("--ntarg"), type="integer", help="number of targets"),
    make_option(c("--nsubj"), type="integer", help="number of subjects"),
    #make_option(c("--pvarHi"), type="integer", help="variance of primes in 'high' condition"),
    #make_option(c("--pvarLo"), type="integer", help="variance of primes in 'low' condition")
    make_option(c("--varianceLst"), type="character", help="comma separated values"),
    make_option(c("--dateStr"), type="character", help="string representing the date")
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
#ncores <<- 2 * n.iter - 1 # This should be one less than is actually 
# available. Reserve one core for the main process.

# These are only needed if generating the data in files ahead of time.

resultDir <- 'Results'
dataDir <- 'GeneratedData'
timingFile <- 'Timing.txt'; timingFileLock <- paste0(timingFile, '.lock')
timingDir <- 'TimingInfo'
# scratchDirHi <- './scratch/HighVarData'
#if (!dir.exists(resultDir)) {dir.create(resultDir)}
#if (!dir.exists(dataDir)) {dir.create(dataDir)}
if (!file.exists(timingFile)) {file.create(timingFile)}
if (!file.exists(timingFileLock)) {file.create(timingFileLock)}
#if (!dir.exists(timingDir)) {dir.create(timingDir)}

nsubj <<- 15
nprim <<- 2
npcat <<- 2
ntarg <<- 2
ntcat <<- 2
nreps <<- 2 # should be even number

svar <<- 1
# pvarLo <<- 1
# pvarHi <<- pvarLo * 2
tvar <<- 1
evar <<- 1
varianceLst <<- "1"
iter.start <<- 1

# Substitute Provided Arguments for Default Values ====

# Overwrite the defaults when they are provided.
for (var in names(args)) {
  assign(var, args[var][[1]])
}

# process the variance list parameter
varianceLst <<- as.numeric(unlist(strsplit(varianceLst, ",")))

# Print Important Parameters for the Logs ====
varLst <- c(
  #"ncores",
  "n.iter",
  "nsubj",
  "nprim",
  "npcat",
  "ntarg",
  "ntcat",
  "nreps",
  "svar",
  #"pvarLo",
  #"pvarHi",
  "tvar",
  "evar",
  "varianceLst",
  "iter.start"
)
for (var in varLst) {
  comm.print(paste0(var, ": ", get(var)))
}
comm.print(paste('Scratch directory:', dateStr))

# Random Seed & Data Directory ====

# Random number considerations. This will make the result reproducible 
# and also ensure that each iteration is reasonably independent.
# RNGkind("L'Ecuyer-CMRG")
comm.set.seed(593065038)
if (!dir.exists(dateStr)) {dir.create(dateStr)}

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
  
  snum <- 1:nsubj
  pnum <- 1:nprim
  tnum <- 1:ntarg
  rnum <- 1:nreps
  
  prej <- rnorm(snum,0,svar)
  basert <- rnorm(snum,0,1)
  subj <- data.frame(snum,prej,basert)
  
  
  # prime differences ====
  pprot <- rep(rnorm(nprim,0,pvar),npcat)
  pcat <- rep(rnorm(npcat,0,1),each=nprim)
  prime <- data.frame(pnum,pcat,pprot)
  
  
  # target differences ====
  tvaln <- rep(rnorm(ntarg,0,tvar),ntcat)
  tcat <- rep(rnorm(ntcat,0,1),each=ntarg)
  target <- data.frame(tnum,tcat,tvaln)
  
  
  # build integrated data file ====
  d <- expand.grid(rnum=rnum,tnum=tnum,pnum=pnum,snum=snum)
  d <- merge(d,target)
  d <- merge(d,prime,by.x="pnum",by.y="pnum")
  d <- merge(d,subj,by.x="snum",by.y="snum")
  
  d$error <- rnorm(nrow(d),0,evar)
  # d$rt <- 600 + 1*(d$pcat*d$tcat) + 1*(d$pcat*d$tcat*d$prej) + 
  #   1*(d$pcat*d$tcat*d$prej*d$tvaln) + 1*(d$pcat*d$tcat*d$prej*d$pprot) + 
  #   1*(d$pcat*d$tcat*d$prej*d$pprot*d$tvaln) + 5*d$error
  
  # This formula was the first attempt at using real-world values.
  # d$rt <- 6.4 + 0.07*(d$prej) + 0.00*(d$pcat) + 0.02*(d$tcat) + 0.02*(d$pprot) + 
  #   0.02*(d$tvaln) + 0.00*(d$prej*d$pcat) + 0.00*(d$prej*d$tcat) + 
  #   0.04*(d$pcat*d$tcat) + 0.04*(d$prej*d$pprot) + 0.02*(d$pprot*d$tcat) + 
  #   0.00*(d$prej*d$tvaln) + 0.00*(d$pcat*d$tvaln) + 0.00*(d$pprot*d$tvaln) + 
  #   0.08*(d$prej*d$pcat*d$tcat) + 0.00*(d$prej*d$pcat*d$tcat*d$tvaln) + 
  #   0.06*(d$prej*d$pcat*d$tcat*d$pprot) +
  #   0.05*(d$pcat*d$tcat*d$prej*d$pprot*d$tvaln) + 0.25*d$error
  
  # These are the real world values based on the heirarchical ordering approach.
  d <- within(d, {
    rt <- 6.4 + 
      .09449*(prej) + .09449*(pcat) + .09449*(tcat) + 
      .09449*(pprot) + .09449*(tvaln) + 
      # two ways
      .08452*(prej*pcat) + 
      .08452*(prej*tcat) + .08452*(prej*pprot) + .08452*(prej*tvaln) + 
      .08452*(pcat*tcat) + .08452*(pcat*pprot) + .08452*(pcat*tvaln) + 
      .08452*(tcat*pprot) + .08452*(tcat*tvaln) + 
      .08452*(pprot*tvaln) + 
      # three ways
      0.07319*(prej*pcat*tcat) + 
      0.07319*(prej*pcat*pprot) + 0.07319*(prej*pcat*tvaln) + 
      0.07319*(prej*tcat*pprot) + 0.07319*(prej*tcat*tvaln) + 
      0.07319*(prej*pprot*tvaln) + 0.07319*(pcat*tcat*pprot) + 
      0.07319*(pcat*tcat*tvaln) + 0.07319*(pcat*pprot*tvaln) + 
      0.07319*(tcat*pprot*tvaln) + 
      # four ways
      0.05976*(prej*pcat*tcat*pprot) + 
      0.05976*(prej*pcat*tcat*tvaln) + 
      0.05976*(prej*pcat*pprot*tvaln) + 
      0.05976*(prej*tcat*pprot*tvaln) + 
      0.05976*(pcat*tcat*pprot*tvaln) + 
      # five way
      .04226*(pcat*tcat*prej*pprot*tvaln) + 
      0.10351*error
    
    # variance decomposition ====
    snum <- as.factor(snum)
    pcat <- as.factor(pcat)
    pnum <- as.factor(pnum)
    tcat <- as.factor(tcat)
    tnum <- as.factor(tnum)
  })
  
  contrasts(d$snum) <- contr.poly
  contrasts(d$pcat) <- contr.poly
  contrasts(d$pnum) <- contr.poly
  contrasts(d$tcat) <- contr.poly
  contrasts(d$tnum) <- contr.poly
  
  return(d)
}

modelFn <- function (d, i=-1) 
{
  initTime <- proc.time()[3]
  m1 <- lm(rt ~ snum * pcat * tcat * pnum * tnum, data=d)
  lmTm <- round((proc.time()[3] - initTime) / 60, 2)
  #print(paste0('Iteration ', i, ' lm() call finished: ', tm, ' minutes'))
  
  initTime <- proc.time()[3]
  my.anova1 <- Anova(m1, type="III", singular.ok = TRUE)
  anovaTm <- round((proc.time()[3] - initTime) / 60, 2)
  #print(paste0('Iteration ', i, ' Anova() call finished: ', tm, ' minutes'))
  numPar <- nrow(my.anova1)
  est <- my.anova1[1:numPar,'Sum Sq'] / my.anova1[1:numPar,'Df']
  names(est) <- rownames(my.anova1)
  
  # split-half reliability
  
  # Sample half of the trials for each participant by replication cell.
  d$half <- 1
  for (s in unique(d$snum)) {
    curRows <- row.names(d[d$snum==s,])
    sampledRows <- sample(curRows, size=length(curRows)/2, replace=F)
    d[sampledRows, 'half'] <- 2
  }
  
  d1 <- data.frame(snum=unique(d$snum))

  rts <- tapply(d$rt, INDEX=list(d$snum, d$half, d$pcat, d$tcat), mean, na.rm=T)
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
  rSh <- cor(d1$half1,d1$half2)
  est['r_sh'] <- rSh

  # parallel forms reliability, must recode replications to only 
  # 2 values (odd and even)
  d$rnumx <- as.numeric(d$rnum) %% 2
  d$rnumx[d$rnumx==0] <- 2
  d1 <- cast(d, snum ~ rnumx + pcat + tcat, mean, value="rt")
  d1$rnum1 <- (d1[,2] - d1[,3]) - (d1[,4] - d1[,5])
  d1$rnum2 <- (d1[,6] - d1[,7]) - (d1[,8] - d1[,9])
  
  # Build the return vector
  est['r_pf'] <- cor(d1$rnum1,d1$rnum2)
  
  rm(m1, my.anova1, numPar, d, d1)
  
  # Also return the timing info
  tming <- list(lmTm=lmTm, anovaTm=anovaTm)
  return(list(est=est, tming=tming))
}

iterFn <- function (i, curPvar) {
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
  # Find out the 1000s group:
  folderGroup <- floor(i / 1000)
  curDataDir <- paste0(dateStr, '/', dataDir, folderGroup)
  curResultDir <- paste0(dateStr, '/', resultDir, folderGroup)
  curTimingDir <- paste0(dateStr, '/', timingDir, folderGroup)
  if (!dir.exists(curDataDir)) {dir.create(curDataDir)}
  if (!dir.exists(curResultDir)) {dir.create(curResultDir)}
  if (!dir.exists(curTimingDir)) {dir.create(curTimingDir)}
  # Save the data for posterity.
  dataFlNm <- paste0(curDataDir, '/DataIter', i, '.RData')
  save(d, file=dataFlNm)
  #print(paste0('Iteration ', i, ' data gen time: ', tm))
  
  #initTm <- proc.time()[3]
  curEst <- modelFn(d, i=i)$est
  curEst['nprim'] <- nprim
  curEst['ntarg'] <- ntarg
  curEst['nreps'] <- nreps
  # Save the result.
  estFlNm <- paste0(curResultDir, '/EstIter', i, '.csv')
  write.table(t(c(curEst, curPvar)), 
              file=estFlNm,
              row.names = F,
              col.names = T,
              sep=','
              )
  tm <- proc.time()[3] - initTm
  #comm.print(paste0('Iteration ', i, ' model time: ', tm))
  
  # Save the timing info in a separate file
  tmDf <- data.frame(iteration=i, time=tm)
  write.csv(tmDf,
            file=paste0(curTimingDir, '/TmIter', i, '.csv'),
            row.names=F
            )
  
  
  rm(d, initTm); gc();
  return(i)
}

 
# Parallelize! the modeling part... ====

init()

# Make a guiding list of variances and iteration numbers.
iterVec <- iter.start:(n.iter*length(varianceLst))
guideMat <- data.frame(
  i=iterVec,
  variance=rep(varianceLst, times=n.iter)[iterVec]
)

time.proc <- system.time({
  id <- get.jid(nrow(guideMat))
  estLst <- lapply(id, 
                     function (i) {
                       curVar <- guideMat$variance[i]
                       iterNum <- guideMat$i[i]
                       return(
                         c(iterFn(iterNum, curPvar=curVar), var=curVar)
                       )
                     }
  )
  #estLst <- unlist(allgather(estLst), recursive=F)
})
comm.print(time.proc)
#file.remove(timingFileLock)

endTime <- suppressWarnings(format(Sys.time(), format='%Y-%m-%d_%H.%M.%S'))
comm.print(paste("Run finished", endTime))

finalize()
