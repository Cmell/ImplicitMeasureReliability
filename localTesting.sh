#!/bin/bash

DATE=$(date +%Y-%m-%d_%H.%M.%S)
mpiexec -np 6 Rscript rxx031mpi.R --n.iter=100 --dateStr=$DATE
Rscript rxxFormResults.R --dateStr=$DATE