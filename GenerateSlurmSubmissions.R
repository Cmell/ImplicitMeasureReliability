baseScript <- '
#!/bin/bash

#SBATCH --partition shas
#SBATCH --job-name jobNameSUB
#SBATCH --time 20:00:00
#SBATCH --nodes numNodesSUB
#SBATCH --ntasks-per-node 24
#SBATCH --qos normal
#SBATCH --output ImpRel_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chme2908@colorado.edu

module purge
module load R/3.3.0
module load openmpi/1.10.2

DATE=$(date +%Y-%m-%d_%H.%M.%S)
let ncores=SLURM_NTASKS-1
echo "Number of cores:" $ncores
mpiexec -np $ncores Rscript rxx031mpi.R --n.iter=10000 --dateStr=$DATE --estFlNm=estFlNmSUB --guideFl=guideFlSUB --group=groupSUB
'

numGroups <- 11

for (g in 1:numGroups) {
  # Args:
  jobName <- paste0("ImpRel", g)
  numNodes <- 10
  estFlNm <- paste0('Results/est_group', g, '.RData')
  guideFl <- "FullGuide.RData"
  group <- g
  
  argLst <- c(
    "jobName",
    "numNodes",
    "estFlNm",
    "guideFl",
    "group"
  )
  
  curScript <- baseScript
  if (!dir.exists("SlurmScripts")) {dir.create("SlurmScripts")}
  for (a in argLst) {
    curScript <- sub(
      paste0(a, "SUB"),
      get(a),
      curScript
    )
  }
  write(curScript,
        file = paste0("./SlurmScripts/rxxGroup", g, ".sh")
  )
}