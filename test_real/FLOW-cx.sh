#!/bin/sh
# this script can be only used CX.
#PJM --rsc-list "rscunit=cx"
#PJM --rsc-list "rscgrp=cx-small"
#PJM --rsc-list "node=4,elapse=24:00:00"
#PJM --mpi "proc=32"
#PJM --name "test.sh"
#PJM -j
#PJM -S

############ configure ###############
NODES=4
CORES=5
PROCS=32
cnf=user.conf
log=log.solver.txt
slv=./solver.exe
############ configure ###############

export OMP_NUM_THREADS=${CORES}
export OMP_STACKSIZE=256000
# export PARALLEL=${PROCS}

module load intel/2022.3

rm -f result/*dmp* result/*mon* result/*geo*
mpiexec.hydra ./$slv < ./$cnf > ./$log

