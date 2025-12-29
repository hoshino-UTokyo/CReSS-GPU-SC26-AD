#!/bin/sh
# Job script for MIYABI supercomputer - postprocessing (unite)
#PBS -q debug-g
#PBS -l select=1
#PBS -W group_list=jh250015
#PBS -l walltime=00:30:00
#PBS -N cress-unite
#PBS -j oe

############ configure ###############
cnf=user.conf
log=log.unite.txt
exe=./unite.exe
############ configure ###############

cd ${PBS_O_WORKDIR}

${exe} < ${cnf} > ${log}
