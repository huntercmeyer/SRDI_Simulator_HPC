#!/bin/bash

sbatch --export=ALL,L2=50,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=60, singlenode.sh
sbatch --export=ALL,L2=55,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=55, singlenode.sh
sbatch --export=ALL,L2=60,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=50, singlenode.sh
sbatch --export=ALL,L2=65,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=45, singlenode.sh
sbatch --export=ALL,L2=70,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=40, singlenode.sh
sbatch --export=ALL,L2=75,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=35, singlenode.sh
sbatch --export=ALL,L2=80,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=30, singlenode.sh
sbatch --export=ALL,L2=85,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=25, singlenode.sh
sbatch --export=ALL,L2=90,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=20, singlenode.sh
sbatch --export=ALL,L2=95,LAMBDA=0.0496,P=1.8,W=120.0,SHAPE=4,L1=15, singlenode.sh
