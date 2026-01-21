############ configure ###############
NODES=1
PROCS=1
CORES=72
cnf=user.conf
DATE=`date '+%Y%m%d%H%M%S'`
log=log.solver.${DATE}.txt
slv=./solver.exe
############ configure ###############

export OMP_NUM_THREADS=${CORES}
export OMP_STACKSIZE=256000

rm -f result/*dmp* result/*mon* result/*geo*
mpirun -np ${PROCS} ${slv} < ${cnf} > ${log} 2> log.solver.${DATE}.err
