#!/bin/bash

for DIST in 350
 do
 for LAMBDA in 0.44
 do
     qsub -v DISTARG=$DIST,LAMBDAARG=$LAMBDA pbs_job_OMP.script -A hpc_neutron100 
 done
done
