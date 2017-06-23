if (!require(parallel)) {
  install.packages('parallel')
  library(parallel)
}

print(detectCores())
