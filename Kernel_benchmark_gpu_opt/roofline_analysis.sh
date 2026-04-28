#!/bin/bash
#
# roofline_analysis.sh - Roofline analysis for CReSS GPU kernel benchmarks
#
# Usage:
#   ./roofline_analysis.sh [kernel_dir ...]
#
# If no arguments, runs on a representative set of trial kernels.
# To run all:
#   ./roofline_analysis.sh [0-9]*_*/
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${SCRIPT_DIR}/roofline_results"
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
SUMMARY="${OUTDIR}/roofline_summary_${TIMESTAMP}.csv"
LOG="${OUTDIR}/roofline_${TIMESTAMP}.log"

# ncu metrics for roofline model
NCU_METRICS="\
sm__sass_thread_inst_executed_op_fadd_pred_on.sum,\
sm__sass_thread_inst_executed_op_fmul_pred_on.sum,\
sm__sass_thread_inst_executed_op_ffma_pred_on.sum,\
sm__sass_thread_inst_executed_op_dadd_pred_on.sum,\
sm__sass_thread_inst_executed_op_dmul_pred_on.sum,\
sm__sass_thread_inst_executed_op_dfma_pred_on.sum,\
dram__bytes.sum,\
lts__t_bytes.sum,\
l1tex__t_bytes.sum,\
gpu__time_duration.sum"

# GH200 120GB theoretical peaks (adjust for different GPUs)
PEAK_FP32_GFLOPS=59300    # ~59.3 TFLOPS
PEAK_FP64_GFLOPS=29600    # ~29.6 TFLOPS
PEAK_BW_GBS=4000           # ~4 TB/s HBM3

# Trial kernels: a mix of memory-bound, compute-bound, and boundary kernels
TRIAL_KERNELS=(
    "067_copy3d_s_copy3d"           # Pure memory copy (baseline)
    "091_diver3d_s_diver3d"         # 3D divergence (stencil)
    "094_diverpiv_s_diverpiv"       # Divergence for pressure iteration
    "230_pgrad_subroutine"          # Pressure gradient
    "308_steppi_s_steppi"           # Pressure time integration
    "313_stepuv_subroutine"         # UV momentum time integration
    "315_stepwi_subroutine"         # W momentum time integration
    "015_advs_subroutine"           # Scalar advection
    "012_advbspi_subroutine"        # Advection for pressure iteration
    "046_buoywsi_s_buoywsi"         # Buoyancy for W
    "115_gaussel_subroutine"        # Gaussian elimination (tridiagonal)
    "303_smoo4s_subroutine"         # 4th-order smoothing
    "236_phy2cnt_s_phy2cnt"         # Physical-to-computational coord
    "327_timeflt_s_timeflt"         # Time filter (Asselin)
    "321_swp2nxt_s_swp2nxt"         # Swap to next timestep
)

# ---------------------------------------------------------------
# Functions
# ---------------------------------------------------------------

log_msg() {
    local msg="[$(date '+%H:%M:%S')] $1"
    echo "$msg"
    echo "$msg" >> "$LOG"
}

run_one_benchmark() {
    local dir_path="$1"
    local bench_name
    bench_name=$(basename "$dir_path")

    local exe="${dir_path}/kernel_benchmark"
    if [ ! -x "$exe" ]; then
        log_msg "SKIP ${bench_name}: no executable"
        return
    fi
    if [ ! -f "${dir_path}/benchmark.conf" ]; then
        log_msg "SKIP ${bench_name}: no benchmark.conf"
        return
    fi

    log_msg "Profiling ${bench_name} ..."

    # Profile ALL kernel launches, then aggregate in Python.
    # ncu replays each kernel multiple times for metrics, so we only
    # capture the first iteration after warmup. But since we don't know
    # how many kernels per iteration, we capture generously and
    # deduplicate by kernel name in the Python aggregator.
    #
    # Strategy: skip warmup launches, capture enough for one iteration.
    # Warmup count is line 3 of benchmark.conf.
    local warmup
    warmup=$(sed -n '3p' "${dir_path}/benchmark.conf" 2>/dev/null || echo "2")

    # Count kernel launches per full run to compute kernels-per-iteration.
    # Use a lightweight ncu invocation (one metric, all launches).
    local total_launches
    total_launches=$(cd "$dir_path" && \
        timeout 60 ncu --metrics gpu__time_duration.sum \
            --csv --kernel-name-base demangled \
            ./kernel_benchmark 2>/dev/null \
        | grep -c '^"[0-9]') || true
    total_launches=${total_launches:-0}

    local bench_iter
    bench_iter=$(sed -n '2p' "${dir_path}/benchmark.conf" 2>/dev/null || echo "10")
    local total_iter=$((warmup + bench_iter))

    local max_kpi=20
    if [ "$total_launches" -gt 0 ] && [ "$total_iter" -gt 0 ]; then
        max_kpi=$(( total_launches / total_iter ))
        [ "$max_kpi" -lt 1 ] && max_kpi=1
    fi
    local skip_count=$((warmup * max_kpi))

    log_msg "  launches=${total_launches}, kpi=${max_kpi}, skip=${skip_count}"

    # Write ncu output to a temp file to avoid pipe/variable issues
    local tmpcsv
    tmpcsv=$(mktemp "${OUTDIR}/ncu_tmp_XXXXXX.csv")

    (cd "$dir_path" && \
        timeout 300 ncu --metrics "${NCU_METRICS}" \
            --csv \
            --kernel-name-base demangled \
            --launch-skip "${skip_count}" \
            --launch-count "${max_kpi}" \
            ./kernel_benchmark 2>/dev/null \
        | grep '^"' > "$tmpcsv") || true

    local nlines
    nlines=$(wc -l < "$tmpcsv")
    if [ "$nlines" -le 1 ]; then
        # Fallback: skip=0, capture first iteration
        log_msg "  Retrying with skip=0 ..."
        (cd "$dir_path" && \
            timeout 300 ncu --metrics "${NCU_METRICS}" \
                --csv \
                --kernel-name-base demangled \
                --launch-skip 0 \
                --launch-count $((max_kpi > 1 ? max_kpi : 20)) \
                ./kernel_benchmark 2>/dev/null \
            | grep '^"' > "$tmpcsv") || true
        nlines=$(wc -l < "$tmpcsv")
    fi

    if [ "$nlines" -le 1 ]; then
        log_msg "  WARN: no ncu data"
        rm -f "$tmpcsv"
        return
    fi

    # Parse and append to summary
    python3 << PYEOF
import csv, sys
from collections import defaultdict

bench_name = "${bench_name}"
tmpcsv = "${tmpcsv}"

with open(tmpcsv) as f:
    lines = f.read().strip().split('\n')

header = None
data_lines = []
for line in lines:
    if '"ID"' in line and '"Metric Name"' in line:
        header = line
    elif line.startswith('"') and '"Command line profiler metrics"' in line:
        data_lines.append(line)

if not header or not data_lines:
    print(f"  WARNING: No parseable data for {bench_name}", file=sys.stderr)
    sys.exit(0)

# Group by kernel name - take only the FIRST occurrence of each kernel
kernels = {}
kernel_order = []
reader = csv.DictReader([header] + data_lines)
for row in reader:
    kname = row["Kernel Name"]
    if kname not in kernels:
        kernels[kname] = {}
        kernel_order.append(kname)
    # Only keep the first set of metrics (first iteration)
    mname = row["Metric Name"]
    if mname not in kernels[kname]:
        mval = row["Metric Value"].replace(",", "")
        try:
            kernels[kname][mname] = float(mval)
        except ValueError:
            kernels[kname][mname] = 0.0

with open("${SUMMARY}", "a") as out:
    for kname in kernel_order:
        metrics = kernels[kname]
        fp32_flop = (
            metrics.get("sm__sass_thread_inst_executed_op_fadd_pred_on.sum", 0)
            + metrics.get("sm__sass_thread_inst_executed_op_fmul_pred_on.sum", 0)
            + 2 * metrics.get("sm__sass_thread_inst_executed_op_ffma_pred_on.sum", 0)
        )
        fp64_flop = (
            metrics.get("sm__sass_thread_inst_executed_op_dadd_pred_on.sum", 0)
            + metrics.get("sm__sass_thread_inst_executed_op_dmul_pred_on.sum", 0)
            + 2 * metrics.get("sm__sass_thread_inst_executed_op_dfma_pred_on.sum", 0)
        )
        total_flop = fp32_flop + fp64_flop

        dram_bytes = metrics.get("dram__bytes.sum", 0)
        l2_bytes = metrics.get("lts__t_bytes.sum", 0)
        l1_bytes = metrics.get("l1tex__t_bytes.sum", 0)
        duration_ns = metrics.get("gpu__time_duration.sum", 0)

        ai_dram = total_flop / dram_bytes if dram_bytes > 0 else 0
        ai_l2 = total_flop / l2_bytes if l2_bytes > 0 else 0
        ai_l1 = total_flop / l1_bytes if l1_bytes > 0 else 0

        duration_s = duration_ns / 1e9
        gflops = (total_flop / 1e9) / duration_s if duration_s > 0 else 0
        bw_gbs = (dram_bytes / 1e9) / duration_s if duration_s > 0 else 0

        # Shorten kernel name
        # e.g. "kernel_benchmark_gpu_stepwi_kernel_stepwi_328" -> "stepwi_328"
        # e.g. "kernel_benchmark_copy3d_kernel_copy3d_184" -> "copy3d_184"
        short = kname
        parts = kname.split("_")
        # Find last occurrence of "kernel" to get the meaningful suffix
        last_kernel_idx = -1
        for i in range(len(parts)):
            if parts[i] == "kernel":
                last_kernel_idx = i
        if last_kernel_idx >= 0 and last_kernel_idx + 1 < len(parts):
            short = "_".join(parts[last_kernel_idx + 1:])

        out.write(f"{bench_name},{short},{fp32_flop:.0f},{fp64_flop:.0f},{total_flop:.0f},"
                  f"{dram_bytes:.0f},{l2_bytes:.0f},{l1_bytes:.0f},"
                  f"{duration_ns:.0f},"
                  f"{ai_dram:.4f},{ai_l2:.4f},{ai_l1:.4f},"
                  f"{gflops:.2f},{bw_gbs:.2f}\n")

PYEOF

    rm -f "$tmpcsv"
    log_msg "  done"
}

# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

mkdir -p "$OUTDIR"
: > "$LOG"

log_msg "Roofline analysis started"
log_msg "Output: ${SUMMARY}"

# Determine target benchmarks
declare -a targets=()
if [ $# -gt 0 ]; then
    for arg in "$@"; do
        arg="${arg%/}"
        if [[ "$arg" != /* ]]; then
            arg="${SCRIPT_DIR}/${arg}"
        fi
        targets+=("$arg")
    done
else
    for k in "${TRIAL_KERNELS[@]}"; do
        if [ -d "${SCRIPT_DIR}/${k}" ]; then
            targets+=("${SCRIPT_DIR}/${k}")
        else
            log_msg "SKIP ${k}: directory not found"
        fi
    done
fi

log_msg "Targets: ${#targets[@]} benchmarks"

# Write CSV header
echo "benchmark,kernel,fp32_flop,fp64_flop,total_flop,dram_bytes,l2_bytes,l1_bytes,duration_ns,ai_dram,ai_l2,ai_l1,gflops,bw_gbs" > "$SUMMARY"

# GPU info
log_msg "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo unknown)"

for target in "${targets[@]}"; do
    run_one_benchmark "$target"
done

log_msg "Profiling complete"
echo ""

# ---------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------

python3 << PYEOF
import csv

summary_file = "${SUMMARY}"
peak_fp32 = ${PEAK_FP32_GFLOPS}
peak_bw = ${PEAK_BW_GBS}
ridge_point = peak_fp32 / peak_bw

print()
width = 135
print("=" * width)
print("ROOFLINE ANALYSIS SUMMARY")
print(f"GPU: GH200 120GB  |  FP32 Peak: {peak_fp32/1000:.1f} TFLOPS  |  HBM3 BW: {peak_bw} GB/s  |  Ridge Point: {ridge_point:.1f} FLOP/Byte")
print("=" * width)

fmt = "{:<45s} {:>12s} {:>10s} {:>10s} {:>10s} {:>12s} {:>12s} {:>10s}"
print(fmt.format(
    "Benchmark / Kernel", "FLOP", "DRAM(MB)", "Time(us)",
    "AI(DRAM)", "GFLOP/s", "BW(GB/s)", "Bound"))
print("-" * width)

rows = []
with open(summary_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append(row)

if not rows:
    print("  (no data)")
else:
    for row in rows:
        bench = row["benchmark"]
        # Shorten: remove leading number prefix
        if "_" in bench:
            bench = bench.split("_", 1)[1]
        kernel = row["kernel"]
        total_flop = float(row["total_flop"])
        dram_bytes = float(row["dram_bytes"])
        duration_ns = float(row["duration_ns"])
        ai_dram = float(row["ai_dram"])
        gflops = float(row["gflops"])
        bw_gbs = float(row["bw_gbs"])

        bound = "MEM" if ai_dram < ridge_point else "COMP"

        if bound == "MEM":
            eff = bw_gbs / peak_bw * 100
            eff_str = f"MEM {eff:.1f}%"
        else:
            eff = gflops / peak_fp32 * 100
            eff_str = f"COMP {eff:.1f}%"

        label = f"{bench} / {kernel}"
        if len(label) > 45:
            label = label[:42] + "..."

        print(fmt.format(
            label,
            f"{total_flop:.2e}",
            f"{dram_bytes/1e6:.1f}",
            f"{duration_ns/1e3:.1f}",
            f"{ai_dram:.4f}",
            f"{gflops:.1f}",
            f"{bw_gbs:.1f}",
            eff_str))

print("-" * width)
print(f"AI < {ridge_point:.1f} => Memory-bound (MEM),  AI >= {ridge_point:.1f} => Compute-bound (COMP)")
print(f"Efficiency: %BW = fraction of peak HBM bandwidth,  %FP = fraction of peak FP32 throughput")
print(f"Results: {summary_file}")
print()
PYEOF
