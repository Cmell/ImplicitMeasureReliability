library(optparse)

opts = 
  optionList = list(
    make_option(c("--nprim"), type="integer", help="number of primes"),
    make_option(c("--n.iter"), type="integer", help="number of iterations"),
    make_option(c("--ncores"), type="integer", help="number of cores to use"),
    make_option(c("--ntarg"), type="integer", help="number of targets"),
    make_option(c("--nsubj"), type="integer", help="number of subjects"),
    make_option(c("--pvarHi"), type="integer", help="variance of primes in 'high' condition"),
    make_option(c("--pvarLo"), type="integer", help="variance of primes in 'low' condition"),
    make_option(c("--varLst"), type="character", help="a comma separated string of values")
  )
optParser = OptionParser(option_list = optionList)
args = parse_args(optParser)

print (args)