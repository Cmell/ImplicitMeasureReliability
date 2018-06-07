#!/bin/bash

DATE=$(date +%Y-%m-%d_%H.%M.%S)
GROUP=1
ESTFL="estDir/est${DATE}_${GROUP}.RData"
mpiexec -np 6 Rscript rxx031mpi.R --n.iter=100 --dateStr=$DATE  --estFlNm=$ESTFL --guideFl=GuideMat1.RData --dataDir=GeneratedDataV2 --resultDir=ResultsV2

# Rscript rxxFormResults.R --dateStr=$DATE
