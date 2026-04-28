#!/bin/sh
# Job script for MIYABI - debug queue (faster scheduling)
#PBS -q debug-g
#PBS -l select=1
#PBS -W group_list=jh250015
#PBS -l walltime=00:30:00
#PBS -N cress-dbg
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

export OMP_NUM_THREADS=${CORES}
export OMP_STACKSIZE=256000

rm -f result/*dmp* result/*mon* result/*geo*
mpirun -np ${PROCS} ${slv} < ${cnf} > ${log}
