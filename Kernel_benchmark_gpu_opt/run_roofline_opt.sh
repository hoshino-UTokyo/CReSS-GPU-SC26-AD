#!/bin/bash
# Roofline analysis for optimized GPU kernels
set -e

BASE_DIR="/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload/Kernel_benchmark_gpu_opt"
OUT_DIR="${BASE_DIR}/roofline_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
CSV="${OUT_DIR}/roofline_summary_${TIMESTAMP}.csv"

mkdir -p "${OUT_DIR}"

METRICS="smsp__sass_thread_inst_executed_op_fadd_pred_on.sum"
METRICS="${METRICS},smsp__sass_thread_inst_executed_op_fmul_pred_on.sum"
METRICS="${METRICS},smsp__sass_thread_inst_executed_op_ffma_pred_on.sum"
METRICS="${METRICS},smsp__sass_thread_inst_executed_op_dadd_pred_on.sum"
METRICS="${METRICS},smsp__sass_thread_inst_executed_op_dmul_pred_on.sum"
METRICS="${METRICS},smsp__sass_thread_inst_executed_op_dfma_pred_on.sum"
METRICS="${METRICS},dram__bytes.sum"
METRICS="${METRICS},lts__t_bytes.sum"
METRICS="${METRICS},l1tex__t_bytes.sum"
METRICS="${METRICS},gpu__time_duration.sum"

echo "benchmark,kernel,fp32_flop,fp64_flop,total_flop,dram_bytes,l2_bytes,l1_bytes,duration_ns,ai_dram,ai_l2,ai_l1,gflops,bw_gbs" > "${CSV}"

KERNELS=(
    "003_adjstni_subroutine"
    "097_eddyvis_subroutine"
    "318_strsten_s_strsten"
    "323_termblk_s_termblk"
)

for kernel_dir in "${KERNELS[@]}"; do
    dir="${BASE_DIR}/${kernel_dir}"
    exe="${dir}/kernel_benchmark"

    if [ ! -x "${exe}" ]; then
        echo "SKIP ${kernel_dir}: no executable"
        continue
    fi

    echo "Profiling ${kernel_dir} ..."
    cd "${dir}"

    # Use ncu with --kernel-id to select 3rd invocation (skip warmup)
    ncu --metrics "${METRICS}" \
        --csv --page raw \
        --target-processes all \
        --kernel-id ::regex:.*:2 \
        "${exe}" 2>/dev/null | python3 -c "
import sys, csv

benchmark = '${kernel_dir}'
reader = csv.reader(sys.stdin)
header = None

for row in reader:
    if not row:
        continue
    # Find header row containing our metric names
    if any('dram__bytes.sum' in cell for cell in row):
        header = row
        continue
    if header is None:
        continue
    if len(row) != len(header):
        continue

    # Build column map
    col = {h: i for i, h in enumerate(header)}

    try:
        kname_key = next(k for k in col if 'Kernel Name' in k or 'launch__kernel_name' in k)
        kname = row[col[kname_key]]
    except (StopIteration, KeyError):
        continue

    def get_val(metric):
        for k, idx in col.items():
            if k == metric:
                try:
                    return float(row[idx].replace(',',''))
                except (ValueError, IndexError):
                    return 0.0
        return 0.0

    fadd = get_val('smsp__sass_thread_inst_executed_op_fadd_pred_on.sum')
    fmul = get_val('smsp__sass_thread_inst_executed_op_fmul_pred_on.sum')
    ffma = get_val('smsp__sass_thread_inst_executed_op_ffma_pred_on.sum')
    dadd = get_val('smsp__sass_thread_inst_executed_op_dadd_pred_on.sum')
    dmul = get_val('smsp__sass_thread_inst_executed_op_dmul_pred_on.sum')
    dfma = get_val('smsp__sass_thread_inst_executed_op_dfma_pred_on.sum')
    dram = get_val('dram__bytes.sum')
    l2 = get_val('lts__t_bytes.sum')
    l1 = get_val('l1tex__t_bytes.sum')
    dur = get_val('gpu__time_duration.sum')

    fp32 = fadd + fmul + 2*ffma
    fp64 = dadd + dmul + 2*dfma
    total = fp32 + fp64

    ai_dram = total/dram if dram > 0 else 0
    ai_l2 = total/l2 if l2 > 0 else 0
    ai_l1 = total/l1 if l1 > 0 else 0
    gflops = total/dur if dur > 0 else 0
    bw_gbs = dram/dur if dur > 0 else 0

    print(f'{benchmark},{kname},{fp32:.0f},{fp64:.0f},{total:.0f},{dram:.0f},{l2:.0f},{l1:.0f},{dur:.0f},{ai_dram:.4f},{ai_l2:.4f},{ai_l1:.4f},{gflops:.2f},{bw_gbs:.2f}')
" >> "${CSV}"

    echo "  done"
done

echo ""
echo "=== Results ==="
cat "${CSV}"
echo ""
echo "Saved to: ${CSV}"
