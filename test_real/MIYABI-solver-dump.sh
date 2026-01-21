#!/bin/sh
# Job script for MIYABI supercomputer - solver with kernel dump
#PBS -q short-g
#PBS -l select=1
#PBS -W group_list=jh250015
#PBS -l walltime=05:00:00
#PBS -N cress-dump
#PBS -j oe

############ configure ###############
NODES=1
PROCS=1
CORES=72
cnf=user.conf
log=log.solver-dump.txt
slv=./solver.exe
############ configure ###############

cd ${PBS_O_WORKDIR}

export OMP_NUM_THREADS=${CORES}
export OMP_STACKSIZE=256000

# Set dump targets for missing kernels (15 total)
export DUMP_TARGETS="bc4news,bcycle,getexner,outdmp,outmxn,phvbcs,phvbcuvw,phvs,phvuvw,setcst2d,siadjst,swadjst,swp2nxt,totals,steps"

rm -f result/*dmp* result/*mon* result/*geo*
mpirun -np ${PROCS} ${slv} < ${cnf} > ${log}
