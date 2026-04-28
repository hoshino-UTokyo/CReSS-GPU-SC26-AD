#!/bin/sh
# Job script for MIYABI supercomputer - ALL CPU test (GPU kernels disabled)
# Purpose: Verify that CPU/OpenMP execution with managed memory produces correct results
#PBS -q short-g
#PBS -l select=1
#PBS -W group_list=jh250015
#PBS -l walltime=05:00:00
#PBS -N cress-allcpu-test
#PBS -j oe

############ configure ###############
NODES=1
PROCS=1
CORES=72
cnf=user.conf
log=log.solver.allcpu_test.txt
slv=./solver.exe
############ configure ###############

cd ${PBS_O_WORKDIR}

export OMP_NUM_THREADS=${CORES}
export OMP_STACKSIZE=256000

rm -f result/*dmp* result/*mon* result/*geo*
mpirun -np ${PROCS} ${slv} < ${cnf} > ${log}
