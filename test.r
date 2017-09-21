args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  n.iter <- as.numeric(args[1])
} else {
  # default value
  n.iter <<- 2 
}

print (n.iter)