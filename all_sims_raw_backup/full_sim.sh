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
cd T0.3
make
md.in = md-00.in
srun -n 1 ./a.out
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd .. 
cp -r T0.3/xvconf T0.6/xvconf
cd T0.6
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
cp -r T0.6/xvconf T0.9/xvconf
cd T0.9
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
cp -r T0.9/xvconf T1.2/xvconf
cd T1.2
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
cp -r T1.2/xvconf T1.5/xvconf
cd T1.5
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
cp -r T1.5/xvconf T1.8/xvconf
cd T1.8
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
cp -r T1.8/xvconf T2.1/xvconf
cd T2.1
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
cp -r T2.1/xvconf T2.4/xvconf
cd T2.4
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
cp -r T2.4/xvconf T2.7/xvconf
cd T2.7
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
cp -r T2.7/xvconf T3.0/xvconf
cd T3.0
make
md.in = md-010203.in
srun -n 1 ./a.out
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-0405.in
srun -n 1 ./a.out
srun -n 1 ./a.out
md.in = md-06.in
srun -n 1 ./a.out
cd ..
