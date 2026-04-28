#!/bin/sh
#PBS -q short-g
#PBS -l select=1
#PBS -W group_list=jh250015
#PBS -l walltime=05:00:00
#PBS -N roofline-all
#PBS -j oe

cd ${PBS_O_WORKDIR}

./roofline_analysis.sh [0-9]*_*/ > roofline_job.log 2>&1
