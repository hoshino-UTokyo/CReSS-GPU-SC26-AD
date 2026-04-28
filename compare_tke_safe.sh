#!/bin/bash
# Compare TKE values - safe_only vs CPU reference
CPU_LOG="test_real/log.solver.20260226223516.txt"
GPU_LOG="test_real/log.solver.txt"

paste <(grep "tkemax" "$CPU_LOG" | awk '{print NR, $3}') \
      <(grep -a "tkemax" "$GPU_LOG" | awk '{print $3}') | \
awk '{
  step=$1; cpu=$2; gpu=$3;
  if(gpu == "") {next}
  if(cpu+0 != 0) {
    rel=((gpu-cpu)/cpu)*100;
    printf "Step %3d: CPU=%s  GPU=%s  rel=%+.4f%%\n", step, cpu, gpu, rel
  } else {
    printf "Step %3d: CPU=%s  GPU=%s\n", step, cpu, gpu
  }
}'
