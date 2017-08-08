#!/bin/bash

#SBATCH --partition shas
#SBATCH --job-name test-job
#SBATCH --time 00:00:05
#SBATCH --nodes 1
#SBATCH --qos debug
#SBATCH --output test-job.out

ml R

Rscript test.r
