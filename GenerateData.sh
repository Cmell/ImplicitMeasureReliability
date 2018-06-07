#!/bin/bash

#SBATCH --partition shas
#SBATCH --job-name GenerateData
#SBATCH --time 03:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --qos normal
#SBATCH --output DataGeneration_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chme2908@colorado.edu

module purge
module load R/3.3.0

Rscript GenerateData.R --n.iter=10000 --nsubj=30 --ntarg=8 --nprim=8 --varianceLst=0,1,2,3,4 --guideFlNm=FullGuide.RData --numGroups=11
