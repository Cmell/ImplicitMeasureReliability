#!/bin/bash

#SBATCH --partition shas
#SBATCH --job-name 36-Iteration-Run
#SBATCH --time 05:00:00
#SBATCH --nodes 3
#SBATCH --ntasks-per-node 24
#SBATCH --qos normal
#SBATCH --output test-job_36Iter20Subj_%j.out

ml R

Rscript rxx031a.r --n.iter=36 --ncores=72  --nsubj=20 --ntarg=8 --nprim=8
