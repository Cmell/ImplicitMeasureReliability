#!/bin/bash

DATE=$(date +%Y-%m-%d_%H.%M.%S)
Rscript GenerateData.R --n.iter=1000 --numGroups=6 --dateStr=$DATE --varianceLst=1,2 --guideFlNm=GuideMat1.RData --dataDir=GeneratedDataV2
