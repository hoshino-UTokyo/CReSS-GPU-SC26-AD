#!/bin/sh
#PBS -q debug-g
#PBS -l select=1
#PBS -W group_list=jh250015
#PBS -l walltime=00:30:00


cd ${PBS_O_WORKDIR}

echo "hello world"
