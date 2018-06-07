pkgLst <- c(
  'optparse',
  'reshape'
)
for (p in pkgLst) {
  library(p, character.only = T, quietly = T)
}

# Date String
suppressWarnings({dateStr <- format(Sys.time(), format='%Y-%m-%d_%H.%M.%S')})

# Working Directory ====
cmDir <- '~chrismellinger/GoogleDrive/ImplicitMeasureReliability/'
corcDir <- '/scratch/summit/chme2908/ImplicitMeasureReliability/'
if (dir.exists(cmDir)) {
  setwd(cmDir)
} else if (dir.exists(corcDir)) {
  setwd(corcDir)
}
print(paste0('Working directory: ', getwd()))

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
    # make_option(c("--iter.start"), type="integer", 
    #             help="
    #             Iteration number to begin at. This is with reference to the
    #             total number of iterations n.iter * length(varianceLst).
    #             "),
    make_option(c("--ncores"), type="integer", help="number of cores to use"),
    make_option(c("--ntarg"), type="integer", help="number of targets"),
    make_option(c("--nsubj"), type="integer", help="number of subjects"),
    make_option(c("--nreps"), type="integer", help="number of subjects"),
    make_option(c("--varianceLst"), type="character", help="comma separated values"),
    make_option(c("--dateStr"), type="character", help="string representing the date"),
    make_option(c("--numGroups"), type="integer", 
                help="number of running groups to include in the guide frame"),
    make_option(c("--guideFlNm"), type="character", 
                help="Filename to save guide file in"),
    make_option(c("--dataDir"), type="character", 
                help="Directory to save generated data in")
    )
optParser = OptionParser(option_list = optionList)
args = parse_args(optParser)

# Check the two required arguments
if (is.null(args$n.iter)) {
  stop(
    paste("Must specify n.iter!")
  )
}
if (is.null(args$guideFlNm)) {
  stop(
    paste("Must specify guideFlNm!")
  )
}
n.iter <- args$n.iter
guideFlNm <- args$guideFlNm

# Default simulation parameters ====

# These are only needed if generating the data in files ahead of time.

dataDir <- 'GeneratedData'

if (!dir.exists(dataDir)) {dir.create(dataDir)}

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
  "tvar",
  "evar",
  "varianceLst",
  "numGroups",
  "guideFlNm",
  "dataDir"
)
for (var in varLst) {
  if (exists(var)) {
    print(paste0(var, ": ", get(var)))
  }
}

# Random Seed & Data Directory ====

# Random number considerations. This will make the result reproducible 
# and also ensure that each iteration is reasonably independent.
set.seed(593065038)

# Function Defintion ====

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
  d$rt <- 6.4 + 0.07*(d$prej) + 0.00*(d$pcat) + 0.02*(d$tcat) + 0.02*(d$pprot) +
    0.02*(d$tvaln) + 0.00*(d$prej*d$pcat) + 0.00*(d$prej*d$tcat) +
    0.04*(d$pcat*d$tcat) + 0.04*(d$prej*d$pprot) + 0.02*(d$pprot*d$tcat) +
    0.00*(d$prej*d$tvaln) + 0.00*(d$pcat*d$tvaln) + 0.00*(d$pprot*d$tvaln) +
    0.08*(d$prej*d$pcat*d$tcat) + 0.00*(d$prej*d$pcat*d$tcat*d$tvaln) +
    0.06*(d$prej*d$pcat*d$tcat*d$pprot) +
    0.05*(d$pcat*d$tcat*d$prej*d$pprot*d$tvaln) + 0.25*d$error
  
  # These are the real world values based on the heirarchical ordering approach.
  d <- within(d, {
    #   rt <- 6.4 + 
    #     .09449*(prej) + .09449*(pcat) + .09449*(tcat) + 
    #     .09449*(pprot) + .09449*(tvaln) + 
    #     # two ways
    #     .08452*(prej*pcat) + 
    #     .08452*(prej*tcat) + .08452*(prej*pprot) + .08452*(prej*tvaln) + 
    #     .08452*(pcat*tcat) + .08452*(pcat*pprot) + .08452*(pcat*tvaln) + 
    #     .08452*(tcat*pprot) + .08452*(tcat*tvaln) + 
    #     .08452*(pprot*tvaln) + 
    #     # three ways
    #     0.07319*(prej*pcat*tcat) + 
    #     0.07319*(prej*pcat*pprot) + 0.07319*(prej*pcat*tvaln) + 
    #     0.07319*(prej*tcat*pprot) + 0.07319*(prej*tcat*tvaln) + 
    #     0.07319*(prej*pprot*tvaln) + 0.07319*(pcat*tcat*pprot) + 
    #     0.07319*(pcat*tcat*tvaln) + 0.07319*(pcat*pprot*tvaln) + 
    #     0.07319*(tcat*pprot*tvaln) + 
    #     # four ways
    #     0.05976*(prej*pcat*tcat*pprot) + 
    #     0.05976*(prej*pcat*tcat*tvaln) + 
    #     0.05976*(prej*pcat*pprot*tvaln) + 
    #     0.05976*(prej*tcat*pprot*tvaln) + 
    #     0.05976*(pcat*tcat*pprot*tvaln) + 
    #     # five way
    #     .04226*(pcat*tcat*prej*pprot*tvaln) + 
    #     0.10351*error
    
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

#  ====

# guide data frame

iterVec <- iter.start:(n.iter*length(varianceLst))
thousandsVec <- (floor((iter.start:n.iter - 1) / 1000)) * 1000
guideMat <- data.frame(
  i=iterVec,
  variance=rep(varianceLst, times=length(iter.start:n.iter))[iterVec],
  thouGroup=rep(thousandsVec, each=length(varianceLst))
)

# Make a filename section that is basically a combination of characteristics
# of the iteration.
guideMat$flNm <- sapply(1:nrow(guideMat), function (r) {
  curVar <- guideMat$variance[r]
  curI <- guideMat$i[r]
  curThou <- guideMat$thouGroup[r]
  flNm <- paste0(dataDir, '/var', curVar, '/set', curThou, '/data', curI, '.RData')
  return(flNm)
})

# Add the groups. At the moment, just divide the rows evenly since each row is 
# run independently.
if (exists("numGroups")) {
  rep(1:numGroups, length.out=nrow(guideMat))
  guideMat$group <- rep(1:numGroups, length.out=nrow(guideMat))
} else {
  guideMat$group <- 1
}

dirLst <- unique(dirname(guideMat$flNm))
print(dirLst)

for (p in dirLst) {
  if (!dir.exists(p)) {dir.create(p, recursive = T)}
}

system.time({
  for (r in 1:nrow(guideMat)) {
    curPvar <- guideMat$variance[r]
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
    save(d, file = guideMat$flNm[r])
  }
})

varLst <- c(varLst, 'guideMat')
varLst <- varLst[sapply(varLst, exists)]
save(list = varLst,  file = guideFlNm)