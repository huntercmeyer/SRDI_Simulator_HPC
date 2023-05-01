#!/bin/bash
#SBATCH -N 1 
#SBATCH -t 10:00:00  
#SBATCH -p workq 
#SBATCH -n 48
#SBATCH -J DG_S8
#SBATCH -A loni_neutronmpg



cd /work/ihidrovo/src/PhaseObjectDualGratings

export OMP_NUM_THREADS=48; ./PhaseStep ${DISTARG} ${LAMBDAARG} ${PHASESTEPARG} 


wait 

exit 

#End-of-file (EOF)
