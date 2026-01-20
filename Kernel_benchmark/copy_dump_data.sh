#!/bin/bash
# Script to copy dump data from test_real/kernel_dump to Kernel_benchmark directories

BASE_DIR="/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload"
DUMP_DIR="$BASE_DIR/test_real/kernel_dump"
BENCH_DIR="$BASE_DIR/Kernel_benchmark"

echo "=== Copying dump data to Kernel_benchmark ==="

for kernel in $(ls "$DUMP_DIR" 2>/dev/null); do
    # Find matching benchmark directory
    dest=$(find "$BENCH_DIR" -maxdepth 1 -type d -name "*_${kernel}_*" -o -name "*_${kernel}" | head -1)

    if [ -z "$dest" ]; then
        echo "WARNING: No matching directory for kernel: $kernel"
        continue
    fi

    if [ -d "$dest/data" ]; then
        # echo "SKIP: $kernel (data already exists)"
        echo "COPY: $kernel -> $(basename $dest)"
        cp -r "$DUMP_DIR/$kernel/"* "$dest/data/"
    else
        echo "COPY: $kernel -> $(basename $dest)"
        mkdir -p "$dest/data"
        cp -r "$DUMP_DIR/$kernel/"* "$dest/data/"
    fi
done

echo ""
echo "=== Summary ==="
echo "Kernels with dump data: $(ls "$DUMP_DIR" | wc -l)"
echo "Kernels with data in Kernel_benchmark: $(find "$BENCH_DIR" -maxdepth 2 -type d -name data | wc -l)"
