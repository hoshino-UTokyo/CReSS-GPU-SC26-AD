#!/bin/bash
cd /work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload/Kernel_benchmark

echo "=== Completed kernel benchmarks ==="
ls -d */kernel_benchmark.f90 2>/dev/null | wc -l

echo ""
echo "=== Incomplete kernels (from omp_profile) ==="

# List of profiled subroutine names
names="strsten turbuvw vspuvw coriuv buoywb rstuvwc sheartke eddyvisj turbtke vsps0 eddydif forcept vsps advbspt smoo4qv vspqv exbcu exbcv diverpih diverpiv rbcw bbcw lbcw exbcss exbcpt exbcq rbcq lbcs rbcs rbcs0 getvdens getexner diagnw termblk setblk nuc1stv nuc1stc collect prodctwg distrpg aggregat melting depsit convers nuc2nd freezing shedding more0q newblk fallblk upwqp upwnp siadjst swadjst phvuvw phvs phvbcuvw phvbcs swp2nxt timeflt kh8uv defomssq eddyvis bc8u bc8v bcten"

for name in $names; do
  for d in *${name}*; do
    if [ -d "$d" ] && [ ! -f "$d/kernel_benchmark.f90" ]; then
      echo "TODO: $d"
    fi
  done
done 2>/dev/null
