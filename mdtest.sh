#!/bin/bash
## Number of tasks you have
#SBATCH --ntasks=1
## how long (in minuites) you want to run your job
#SBATCH --time=60
## output file 

## your job name
#SBATCH --job-name=T0.3             
## change directory
#SBATCH --chdir=.
cd MD2
make
srun -n 1 ./a.out
cd .. 
cp -r MD2/xvconf MD3/xvconf
cd MD3
make
srun -n 1 ./a.out
