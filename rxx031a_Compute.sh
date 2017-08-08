#!/bin/bash

#SBATCH --partition shas
#SBATCH --job-name 32-Iteration-Run
#SBATCH --time 01:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 24
#SBATCH --qos normal
#SBATCH --output 32-job_%j.out

ml R

Rscript rxx031a.r
