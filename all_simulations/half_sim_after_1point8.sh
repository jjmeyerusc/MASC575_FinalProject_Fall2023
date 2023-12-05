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
cp -r T1.8/xvconf T2.1/xvconf
cd T2.1
make
cp -r md-010203.in md.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
cp -r md-0405.in md.in
srun -n 1 ./a.out
srun -n 1 ./a.out
cp -r md-06.in md.in
srun -n 1 ./a.out
cd ..
cp -r T2.1/xvconf T2.4/xvconf
cd T2.4
make
cp -r md-010203.in md.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
cp -r md-0405.in md.in
srun -n 1 ./a.out
srun -n 1 ./a.out
cp -r md-06.in md.in
srun -n 1 ./a.out
cd ..
cp -r T2.4/xvconf T2.7/xvconf
cd T2.7
make
cp -r md-010203.in md.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
cp -r md-0405.in md.in
srun -n 1 ./a.out
srun -n 1 ./a.out
cp -r md-06.in md.in
srun -n 1 ./a.out
cd ..
cp -r T2.7/xvconf T3.0/xvconf
cd T3.0
make
cp -r md-010203.in md.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
cp -r md-0405.in md.in
srun -n 1 ./a.out
srun -n 1 ./a.out
cp -r md-06.in md.in
srun -n 1 ./a.out
cd ..
