#!/bin/bash
#SBATCH -N 1 
#SBATCH -t 10:00:00  
#SBATCH -p workq
#SBATCH -n 48
#SBATCH -J MPGWithObject
#SBATCH -A hpc_mpgxray



#cd /work/hmeyer5/src/PhaseObjectSingleMPG

#export OMP_NUM_THREADS=48; ./RunGrating 250 .44 1.0 300.0 4 50
export OMP_NUM_THREADS=48; ./RunGrating ${L2} ${LAMBDA} ${P} ${W} ${SHAPE} ${L1}


wait 

exit 

#End-of-file (EOF)
