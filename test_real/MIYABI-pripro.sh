#!/bin/sh
# Job script for MIYABI supercomputer - preprocessing (gridata)
#PBS -q debug-g
#PBS -l select=1
#PBS -W group_list=jh250015
#PBS -l walltime=00:30:00
#PBS -N cress-pripro
#PBS -j oe

############ configure ###############
cnf=user.conf
log=log.gridata.txt
exe=./gridata.exe
############ configure ###############

cd ${PBS_O_WORKDIR}

rm -f result/*gpv*
${exe} < ${cnf} > ${log}
