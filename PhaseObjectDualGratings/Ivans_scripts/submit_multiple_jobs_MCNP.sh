#!/bin/bash

for DIST in 12
 do
 for LAMBDA in 0.44
 do
     sbatch -v DISTARG=$DIST,LAMBDAARG=$LAMBDA singlenode.sh 
 done
done
