# Get Args ====
library(optparse)
opts = 
  optionList = list(
    make_option(c("--dateStr"), type="character", help="string representing the date")
    )
optParser = OptionParser(option_list = optionList)
args = parse_args(optParser)

# Overwrite variables.
for (var in names(args)) {
  assign(var, args[var][[1]])
}

# Functions ====

finishCalcs <- function (est) {
  est <- within(est, {
    var.resid <- ms.resid
    var.snpctc <- 
      (ms.snpctc - ms.snpctn - ms.sntcpn + ms.snpntn) / 
      (nprim*ntarg*nreps)
    var.snpctn <- (ms.snpctn - ms.snpntn) / (nprim*nreps)
    var.sntcpn <- (ms.sntcpn - ms.snpntn) / (ntarg*nreps)
    var.snpntn <- (ms.snpntn - ms.resid) / (nreps)
    
    # recode negative variances to zero (this will bias the estimates but is 
    # necessary to avoid negative reliabilities)
    var.snpctc[var.snpctc<0] <- 0
    var.snpctn[var.snpctn<0] <- 0
    var.sntcpn[var.sntcpn<0] <- 0
    var.snpntn[var.snpntn<0] <- 0
    
    rxxmse <- (ms.snpctc-ms.resid) / ms.snpctc
    
    rxxvar <- 
      (var.snpctc + (var.snpctn / ntarg) + (var.sntcpn/ nprim) + 
         (var.snpntn / (nprim*ntarg))
      ) / 
      (var.snpctc + (var.snpctn / ntarg) + (var.sntcpn / nprim) + 
         (var.snpntn / (nprim*ntarg)) + var.resid / (nprim*ntarg*nreps)
      )
    
    rxxvar.prand <- 
      var.snpctc / 
      (var.snpctc + (var.snpctn / ntarg) + (var.sntcpn / nprim) + 
         (var.snpntn / (nprim*ntarg)) + var.resid / (nprim*ntarg*nreps)
      )
  })
  
  return(est)
}

renameEstCols <- function (est) {
  colnames(est) <- c(
    "ms.int",
    "ms.sn",
    "ms.pc",
    "ms.tc",
    "ms.pn",
    "ms.tn",
    
    "ms.snpc",
    "ms.sntc",
    "ms.pctc",
    "ms.snpn",
    "ms.pcpn",
    "ms.tcpn",
    "ms.sntn",
    "ms.pctn",
    "ms.tctn",
    "ms.pntn",
    
    "ms.snpctc",
    "ms.snpcpn",
    "ms.sntcpn",
    "ms.pctcpn",
    "ms.snpctn",
    "ms.sntctn",
    "ms.pctctn",
    "ms.snpntn",
    "ms.pcpntn",
    "ms.tcpntn",
    
    "ms.snpctcpn",
    "ms.snpctctn",
    "ms.snpcpntn",
    "ms.sntcpntn",
    "ms.pctcpntn",
    
    "ms.snpctcpntn",
    "ms.resid",
    "r_sh","r_pf", "nprim", "ntarg", "nreps", "var"
  )
  return(est)
}

# Working Directory ====

cmDir <- '~chrismellinger/GoogleDrive/ImplicitMeasureReliability/'
corcDir <- '/scratch/summit/chme2908/ImplicitMeasureReliability/'
if (dir.exists(cmDir)) {
  setwd(cmDir)
} else if (dir.exists(corcDir)) {
  setwd(corcDir)
}

resultDir <- 'Results'

# Make the est dataframe. ====
# Load Data
flLst <- list.files(path=dateStr, full.names = T, recursive = T)
flLst <- grep(paste0(resultDir, ".*csv"), flLst, value=T)
est <- read.csv(flLst[[1]], header=T)
for (fl in flLst[2:length(flLst)]) {
  curEst <- read.csv(fl, header=T)
  est <- rbind(est, curEst)
}

est <- renameEstCols(est)

# Perform calculations. ====

est <- finishCalcs(est)

# Save all the things. ====

print(paste('Results formed for', dateStr))

save(est,
     file=paste0('est_', dateStr,'.RData')
)