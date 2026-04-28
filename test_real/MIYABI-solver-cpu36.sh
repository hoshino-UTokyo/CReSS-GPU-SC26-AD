#!/bin/sh
# Job script for MIYABI - CPU-only solver, 36 timesteps
#PBS -q short-g
#PBS -l select=1
#PBS -W group_list=jh250015
#PBS -l walltime=05:00:00
#PBS -N cress-cpu36
#PBS -j oe

############ configure ###############
NODES=1
PROCS=1
CORES=72
cnf=user_cpu36.conf
log=log.solver_cpu36.txt
############ configure ###############

cd ${PBS_O_WORKDIR}
TOPDIR=$(cd .. && pwd)

# Step 1: Build CPU-only solver
echo "=== Building CPU-only solver ==="
cd ${TOPDIR}/Src

# Ensure Src/Src symlinks exist for Makefile
mkdir -p Src 2>/dev/null
ln -sf ../Make_solver Src/Make_solver 2>/dev/null
ln -sf ../Make_object Src/Make_object 2>/dev/null
ln -sf ../Make_module Src/Make_module 2>/dev/null

# Use CPU-only compile.conf
cp ${TOPDIR}/compile.conf.cpu compile.conf

# Clean and rebuild
make clean TARGET=solver 2>&1 | tail -3
make solver TARGET=solver 2>&1 | tail -20
echo "=== Build complete ==="
ls -la solver.exe

# Step 2: Copy to test_real and run
cd ${PBS_O_WORKDIR}
cp ${TOPDIR}/Src/solver.exe ./solver_cpu.exe

echo "=== Running CPU solver (36 steps) ==="
export OMP_NUM_THREADS=${CORES}
export OMP_STACKSIZE=256000

rm -f result/*dmp* result/*mon* result/*geo*
mpirun -np ${PROCS} ./solver_cpu.exe < ${cnf} > ${log}

echo "=== Done ==="
tail -30 ${log}
