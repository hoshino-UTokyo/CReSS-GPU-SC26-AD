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

# rm -f result/*dmp* result/*mon* result/*geo*
# mpirun -np ${PROCS} ${slv} < ${cnf} > ${log} 2> log.solver.${DATE}.err

# rm -f result/*dmp* result/*mon* result/*geo*
# time mpirun -np ${PROCS} solver_cpu_no_acc.exe < ${cnf} > log.solver_noacc.${DATE}.txt 2> log.solver_noacc.${DATE}.err


export NVCOMPILER_ACC_TIME=1
rm -f result/*dmp* result/*mon* result/*geo*
time mpirun -np ${PROCS} solver.exe < ${cnf} > log.solver_gpu.${DATE}.txt 2> log.solver_gpu.${DATE}.err

