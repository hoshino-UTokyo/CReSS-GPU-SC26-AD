#!/bin/sh
# this script can be only used FX10.
#PJM --rsc-list "rscunit=fx"
#PJM --rsc-list "rscgrp=fx-small"
#PJM --rsc-list "node=4,elapse=24:00:00"
#PJM --mpi "proc=32"
#PJM --name "test.sh"
#PJM -j
#PJM -S

############ configure ###############
NODES=4
CORES=6
PROCS=32
cnf=user.conf
log=log.solver.txt
slv=./solver.exe
############ configure ###############

export OMP_NUM_THREADS=${CORES}
export OMP_STACKSIZE=256000
export PARALLEL=${PROCS}

rm -f result/*dmp* result/*mon* result/*geo*
mpiexec -stdin ./$cnf -of-proc ./$log $slv -Wl,-Lu -Wl,-T

