#!/bin/bash

for DIST in 12
 do
   for LAMBDA in 0.44 
 do
     for PHASESTEP in 0.72 0.96 1.2 1.44 1.68 1.92 2.16 0.0 0.24 0.48 2.4 
      do
      sbatch --export=ALL,DISTARG=$DIST,LAMBDAARG=$LAMBDA,PHASESTEPARG=$PHASESTEP singlenode.sh 
  done
 done
done
