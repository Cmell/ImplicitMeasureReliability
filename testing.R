library(optparse)

opts = 
optionList = list(
  make_option(c("--nprim"), type="integer", help="number of primes"),
  make_option(c("--n.iter"), type="integer", help="number of iterations"),
  make_option(c("--ntarg"), type="integer", help="number of targets"),
  make_option(c("--nsubj"), type="integer", help="number of subjects")
)
optParser = OptionParser(option_list = optionList)
args = parse_args(optParser)
print(names(args))

for (var in names(args)) {
  assign(var, args[var][[1]])
}

for (var in names(args)) {
  print(paste(var, ": ", get(var)))
}