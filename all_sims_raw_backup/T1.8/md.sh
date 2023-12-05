#!/bin/bash
## Number of tasks you have
#SBATCH --ntasks=1
## how long (in minuites) you want to run your job
#SBATCH --time=60
## output file 
#SBATCH --output=md.out
## your job name
#SBATCH --job-name=T0.3             
## change directory
#SBATCH --chdir=.

srun -n 1 ./a.out
