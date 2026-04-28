#!/bin/sh
# Job script for MIYABI supercomputer - solver
#PBS -q short-g
#PBS -l select=1
#PBS -W group_list=jh260061
#PBS -l walltime=05:00:00
#PBS -N cress-solver
#PBS -j oe

############ configure ###############
NODES=1
PROCS=1
CORES=72
cnf=user.conf
log=log.solver.txt
slv=./solver.exe
############ configure ###############

cd ${PBS_O_WORKDIR}

echo "hello"
