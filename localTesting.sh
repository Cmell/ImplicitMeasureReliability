#!/bin/bash

DATE=$(date +%Y-%m-%d_%H.%M.%S)
mpiexec -np 6 Rscript rxx031mpi.R --n.iter=100 --dateStr=$DATE --varianceLst=0,1,2 --nreps=4
Rscript rxxFormResults.R --dateStr=$DATE
