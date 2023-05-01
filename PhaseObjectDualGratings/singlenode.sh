#!/bin/bash
#SBATCH -N 1 
#SBATCH -t 10:00:00  
#SBATCH -p workq 
#SBATCH -n 48
#SBATCH -J DualGrating
#SBATCH -A loni_neutrompg2



cd /work/hmeyer5/src/PhaseObjectDualGratings

#export OMP_NUM_THREADS=48; ./RunDualGratings ${DISTARG} ${LAMBDAARG} ${
# D (mm), Lambda (nm), L1 (cm)
#export OMP_NUM_THREADS=48; ./RunDualGratings 4 0.056 50
export OMP_NUM_THREADS=48; ./RunDualGratings ${D} ${LAMBDA} ${L1} ${PHASESTEPINDEX} ${ACCEPTANCEANGLE}

wait 

exit 

#End-of-file (EOF)
