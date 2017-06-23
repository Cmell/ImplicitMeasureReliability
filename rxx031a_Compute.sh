#!/bin/bash

#SBATCH --partition shas
#SBATCH --job-name test-job
#SBATCH --time 00:10:00
#SBATCH --nodes 1
#SBATCH --qos debug
#SBATCH --output test-job_%j.out

ml R

Rscript rxx031a.r
