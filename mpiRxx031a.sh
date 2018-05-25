#!/bin/bash

#SBATCH --partition shas
#SBATCH --job-name Run_3-MPI
#SBATCH --time 20:00:00
#SBATCH --nodes 90
#SBATCH --ntasks-per-node 24
#SBATCH --qos normal
#SBATCH --output mpiRun_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chme2908@colorado.edu

module purge
module load R/3.3.0
module load openmpi/1.10.2

DATE=$(date +%Y-%m-%d_%H.%M.%S)
let ncores=SLURM_NTASKS-1
echo "Number of cores:" $ncores
mpiexec -np $ncores Rscript rxx031mpi.R --n.iter=10000 --nsubj=30 --ntarg=8 --nprim=8 --varianceLst=0,1,2,3,4 --dateStr=$DATE
# Rscript rxxFormResults.R --dateStr=$DATE
