#!/bin/bash

#SBATCH --job-name test-job
#SBATCH --time 00:00:10
#SBATCH --nodes 6
#SBATCH --output test-job.out

module load slurm/summit
module load gcc
module load intel
module load R

RScript test.r
