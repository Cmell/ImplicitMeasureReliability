#!/bin/bash

#SBATCH --partition shas
#SBATCH --job-name Run_3-MPI
#SBATCH --time 02:00:00
#SBATCH --nodes 2
#SBATCH --ntasks-per-node 24
#SBATCH --qos normal
#SBATCH --output mpiTesting_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chme2908@colorado.edu

module purge
module load R/3.3.0
module load openmpi/1.10.2

DATE=$(date +%Y-%m-%d_%H.%M.%S)
let ncores=SLURM_NTASKS
echo "number of cores:" $ncores
mpiexec -np $ncores Rscript rxx031mpi.R --n.iter=48 --nsubj=30 --ntarg=8 --nprim=8 --varianceLst=1,2 --dateStr=$DATE
# Rscript rxxFormResults.R --dateStr=$DATE
