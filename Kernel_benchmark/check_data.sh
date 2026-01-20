#!/bin/bash
cd /work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload/Kernel_benchmark

echo "=== Kernels with dump data but missing data in Kernel_benchmark ==="
for k in advbspi advs bc4news bcycle buoywsi copy3d diagni disptke diver2d diver3d pgrad phy2cnt steppi steps stepuv stepwi vbcu vbcv vbcw vbcwc; do
  for dest_dir in *${k}*; do
    if [ -d "$dest_dir" ]; then
      if [ ! -d "$dest_dir/data" ]; then
        echo "NEED_COPY: $dest_dir"
      fi
    fi
  done
done
