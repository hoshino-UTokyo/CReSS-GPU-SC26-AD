#!/usr/bin/env python3
"""Roofline comparison: original vs optimized GPU kernels."""

import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# GH200 specs
PEAK_FP32 = 59300
PEAK_FP64 = 29600
PEAK_BW = 4000

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG_CSV = os.path.join(BASE, "../../Kernel_benchmark_gpu/roofline_results/roofline_summary_20260323_233343.csv")

# Find latest opt CSV
opt_csvs = sorted(glob.glob(os.path.join(BASE, "roofline_summary_*.csv")))
OPT_CSV = opt_csvs[-1] if opt_csvs else None

def load_csv(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            ai = float(r['ai_dram'])
            gflops = float(r['gflops'])
            bw = float(r['bw_gbs'])
            rows.append({
                'benchmark': r['benchmark'],
                'kernel': r['kernel'],
                'ai_dram': ai,
                'gflops': gflops,
                'bw_gbs': bw,
                'total_flop': int(float(r['total_flop'])),
                'duration_ns': int(float(r['duration_ns'])),
            })
    return rows

def get_benchmark_time(rows, benchmark):
    """Total time for a benchmark (sum of all kernel durations)."""
    total = sum(r['duration_ns'] for r in rows if r['benchmark'] == benchmark)
    return total

def short_name(benchmark):
    name = benchmark.split('_', 1)[-1] if '_' in benchmark else benchmark
    return name.split('_subroutine')[0].split('_s_')[0]


# ============================================================
# Plot 1: Roofline with optimized kernels highlighted
# ============================================================
def plot_roofline_comparison(orig_data, opt_data):
    fig, ax = plt.subplots(figsize=(14, 9))

    # Roofline line
    ai_range = np.logspace(-3, 3, 500)
    roofline = np.minimum(PEAK_FP64, PEAK_BW * ai_range)
    ax.plot(ai_range, roofline, 'k-', linewidth=2.5,
            label=f'Roofline (FP64={PEAK_FP64} GFLOPS, BW={PEAK_BW} GB/s)')

    # Ridge point
    ridge_ai = PEAK_FP64 / PEAK_BW
    ax.axvline(ridge_ai, color='gray', linestyle=':', alpha=0.5)

    # Original kernels (gray background)
    orig_valid = [d for d in orig_data if d['ai_dram'] > 0 and d['gflops'] > 0]
    ax.scatter([d['ai_dram'] for d in orig_valid],
               [d['gflops'] for d in orig_valid],
               c='lightgray', s=30, alpha=0.4, edgecolors='gray', linewidth=0.3,
               zorder=3, label=f'All kernels (N={len(orig_valid)})')

    # Optimized kernels
    opt_valid = [d for d in opt_data if d['ai_dram'] > 0 and d['gflops'] > 0]
    colors = {'adjstni': 'red', 'eddyvis': 'blue', 'strsten': 'green', 'termblk': 'orange'}
    markers = {'adjstni': 'o', 'eddyvis': 's', 'strsten': 'D', 'termblk': '^'}

    for d in opt_valid:
        sname = short_name(d['benchmark'])
        color = colors.get(sname, 'purple')
        marker = markers.get(sname, 'o')
        ax.scatter(d['ai_dram'], d['gflops'], c=color, s=120, marker=marker,
                   edgecolors='black', linewidth=1, zorder=10)

    # Draw arrows from original to optimized for each benchmark
    opt_benchmarks = set(d['benchmark'] for d in opt_data)
    for bm in opt_benchmarks:
        orig_kernels = [d for d in orig_data if d['benchmark'] == bm and d['ai_dram'] > 0 and d['gflops'] > 0]
        opt_kernels = [d for d in opt_data if d['benchmark'] == bm and d['ai_dram'] > 0 and d['gflops'] > 0]

        for ok in opt_kernels:
            # Find matching original kernel by name
            matches = [o for o in orig_kernels if o['kernel'] and ok['kernel'] and
                       o['kernel'].split('_')[-1] == ok['kernel'].split('_')[-1]]
            if not matches:
                matches = orig_kernels

    # Legend for optimized kernels
    for sname in ['adjstni', 'eddyvis', 'strsten', 'termblk']:
        ax.scatter([], [], c=colors[sname], s=100, marker=markers[sname],
                   edgecolors='black', linewidth=1, label=f'{sname} (optimized)')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Arithmetic Intensity (FLOP/Byte)', fontsize=12)
    ax.set_ylabel('Performance (GFLOPS)', fontsize=12)
    ax.set_title('CReSS GPU Kernel Roofline - Optimized Kernels Highlighted', fontsize=14)
    ax.set_xlim(1e-3, 1e2)
    ax.set_ylim(1e-1, 1e5)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(loc='upper left', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(BASE, 'roofline_comparison.png'), dpi=200, bbox_inches='tight')
    print("Saved: roofline_comparison.png")
    plt.close()


# ============================================================
# Plot 2: Time comparison bar chart
# ============================================================
def plot_time_comparison(orig_data, opt_data):
    # Benchmark times from CSV data
    orig_gpu_dir = os.path.join(BASE, "../../Kernel_benchmark_gpu")
    opt_gpu_dir = os.path.join(BASE, "..")

    benchmarks = ['003_adjstni_subroutine', '097_eddyvis_subroutine',
                  '318_strsten_s_strsten', '323_termblk_s_termblk']

    # Use benchmark execution times from the actual runs
    orig_times = {
        '003_adjstni_subroutine': 1.946,
        '097_eddyvis_subroutine': 1.872,
        '318_strsten_s_strsten': 5.272,
        '323_termblk_s_termblk': 2.721,
    }
    opt_times = {
        '003_adjstni_subroutine': 1.989,
        '097_eddyvis_subroutine': 1.872,  # reverted to original
        '318_strsten_s_strsten': 5.218,
        '323_termblk_s_termblk': 2.735,
    }

    names = [short_name(b) for b in benchmarks]
    orig_vals = [orig_times[b] for b in benchmarks]
    opt_vals = [opt_times[b] for b in benchmarks]
    speedups = [o / p if p > 0 else 1.0 for o, p in zip(orig_vals, opt_vals)]

    x = np.arange(len(names))
    width = 0.35

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Bar chart
    bars1 = ax1.bar(x - width/2, orig_vals, width, label='Original', color='steelblue', edgecolor='white')
    bars2 = ax1.bar(x + width/2, opt_vals, width, label='Optimized', color='coral', edgecolor='white')

    ax1.set_xlabel('Kernel', fontsize=12)
    ax1.set_ylabel('Execution Time (ms)', fontsize=12)
    ax1.set_title('Kernel Execution Time Comparison', fontsize=13)
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, fontsize=10)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, axis='y')

    for bar, val in zip(bars1, orig_vals):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                 f'{val:.3f}', ha='center', va='bottom', fontsize=8)
    for bar, val in zip(bars2, opt_vals):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                 f'{val:.3f}', ha='center', va='bottom', fontsize=8)

    # Speedup chart
    colors_sp = ['green' if s > 1.01 else 'orange' if s > 0.99 else 'red' for s in speedups]
    bars3 = ax2.bar(x, speedups, 0.5, color=colors_sp, edgecolor='white')
    ax2.axhline(1.0, color='black', linestyle='--', linewidth=1)
    ax2.set_xlabel('Kernel', fontsize=12)
    ax2.set_ylabel('Speedup (original/optimized)', fontsize=12)
    ax2.set_title('Optimization Speedup', fontsize=13)
    ax2.set_xticks(x)
    ax2.set_xticklabels(names, fontsize=10)
    ax2.set_ylim(0.9, 1.15)
    ax2.grid(True, alpha=0.3, axis='y')

    for bar, s in zip(bars3, speedups):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                 f'{s:.3f}x', ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(BASE, 'time_comparison.png'), dpi=200, bbox_inches='tight')
    print("Saved: time_comparison.png")
    plt.close()


# ============================================================
# Plot 3: Bandwidth utilization comparison
# ============================================================
def plot_bandwidth_comparison(orig_data, opt_data):
    benchmarks = ['003_adjstni_subroutine', '097_eddyvis_subroutine',
                  '318_strsten_s_strsten', '323_termblk_s_termblk']
    names = [short_name(b) for b in benchmarks]

    fig, ax = plt.subplots(figsize=(12, 7))

    for i, bm in enumerate(benchmarks):
        orig_bws = [d['bw_gbs'] for d in orig_data
                    if d['benchmark'] == bm and d['bw_gbs'] > 0]
        opt_bws = [d['bw_gbs'] for d in opt_data
                   if d['benchmark'] == bm and d['bw_gbs'] > 0]

        if orig_bws:
            ax.scatter([i - 0.15] * len(orig_bws), orig_bws,
                       c='steelblue', s=80, marker='o', edgecolors='black', linewidth=0.5,
                       zorder=5, label='Original' if i == 0 else '')
        if opt_bws:
            ax.scatter([i + 0.15] * len(opt_bws), opt_bws,
                       c='coral', s=80, marker='s', edgecolors='black', linewidth=0.5,
                       zorder=5, label='Optimized' if i == 0 else '')

    ax.axhline(PEAK_BW, color='red', linestyle='--', linewidth=1.5, alpha=0.5,
               label=f'Peak BW ({PEAK_BW} GB/s)')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=11)
    ax.set_ylabel('DRAM Bandwidth (GB/s)', fontsize=12)
    ax.set_title('DRAM Bandwidth Utilization: Original vs Optimized', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, PEAK_BW * 1.1)
    plt.tight_layout()
    plt.savefig(os.path.join(BASE, 'bandwidth_comparison.png'), dpi=200, bbox_inches='tight')
    print("Saved: bandwidth_comparison.png")
    plt.close()


# ============================================================
# Plot 4: Full roofline with optimization summary
# ============================================================
def plot_full_roofline_with_summary(orig_data, opt_data):
    """Full roofline with optimization summary table."""
    fig, (ax, ax_table) = plt.subplots(1, 2, figsize=(18, 9),
                                        gridspec_kw={'width_ratios': [2, 1]})

    # Roofline
    ai_range = np.logspace(-3, 3, 500)
    roofline = np.minimum(PEAK_FP64, PEAK_BW * ai_range)
    ax.plot(ai_range, roofline, 'k-', linewidth=2.5)

    # All original
    orig_valid = [d for d in orig_data if d['ai_dram'] > 0 and d['gflops'] > 0]
    effs = []
    for d in orig_valid:
        roof = min(PEAK_FP64, PEAK_BW * d['ai_dram'])
        effs.append(d['gflops'] / roof if roof > 0 else 0)

    sc = ax.scatter([d['ai_dram'] for d in orig_valid],
                    [d['gflops'] for d in orig_valid],
                    c=effs, cmap='RdYlGn', s=40, alpha=0.6,
                    edgecolors='black', linewidth=0.3, vmin=0, vmax=1, zorder=3)

    # Optimized kernels (larger markers)
    opt_valid = [d for d in opt_data if d['ai_dram'] > 0 and d['gflops'] > 0]
    for d in opt_valid:
        roof = min(PEAK_FP64, PEAK_BW * d['ai_dram'])
        eff = d['gflops'] / roof if roof > 0 else 0
        ax.scatter(d['ai_dram'], d['gflops'], c='blue', s=150,
                   marker='*', edgecolors='black', linewidth=1, zorder=10)
        sname = short_name(d['benchmark'])
        ax.annotate(sname, (d['ai_dram'], d['gflops']),
                    fontsize=7, fontweight='bold', color='blue',
                    xytext=(8, 5), textcoords='offset points')

    cbar = plt.colorbar(sc, ax=ax, shrink=0.6)
    cbar.set_label('Roofline Efficiency', fontsize=10)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Arithmetic Intensity (FLOP/Byte)', fontsize=12)
    ax.set_ylabel('Performance (GFLOPS)', fontsize=12)
    ax.set_title('CReSS GPU Kernel Roofline Model (GH200)', fontsize=14)
    ax.set_xlim(1e-3, 1e2)
    ax.set_ylim(1e-1, 1e5)
    ax.grid(True, alpha=0.3, which='both')

    # Summary table
    ax_table.axis('off')
    ax_table.set_title('Optimization Summary', fontsize=14, pad=20)

    table_data = [
        ['Kernel', 'Orig (ms)', 'Opt (ms)', 'Speedup', 'Optimizations'],
        ['adjstni', '1.946', '1.989', '0.98x', 'collapse(3)+fusion(2→1)'],
        ['strsten', '5.272', '5.218', '1.01x', 'collapse(3)+fusion(8→2)'],
        ['termblk', '2.721', '2.735', '0.99x', 'fusion(2→1)+merge haiopt'],
        ['eddyvis', '1.872', '1.872', '1.00x', 'kappa param fix (bug)'],
        ['', '', '', '', ''],
        ['pgrad*', '4.026', '2.240', '1.80x', 'tmp3 inline+fusion(6→2)'],
        ['gaussel*', '4.674', '2.556', '1.83x', 'k-loop consolidation'],
        ['stepwi*', '3.885', '3.210', '1.21x', 'tmp inline+fusion(5→2)'],
    ]

    table = ax_table.table(cellText=table_data, loc='center',
                           cellLoc='center', colWidths=[0.18, 0.14, 0.14, 0.12, 0.42])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)

    # Header formatting
    for j in range(5):
        table[0, j].set_facecolor('#4472C4')
        table[0, j].set_text_props(color='white', fontweight='bold')

    # Highlight speedups
    for i in range(1, len(table_data)):
        for j in range(5):
            if i <= 4:
                table[i, j].set_facecolor('#E2EFDA')
            elif i == 5:
                table[i, j].set_facecolor('white')
                table[i, j].set_edgecolor('white')
            else:
                table[i, j].set_facecolor('#D6E4F0')

    ax_table.text(0.5, 0.08, '* Previously optimized kernels (for reference)',
                  transform=ax_table.transAxes, ha='center', fontsize=8, style='italic', color='gray')
    ax_table.text(0.5, 0.02,
                  'Note: Memory-bound kernels (AI < ridge point 7.4) are limited by bandwidth.\n'
                  'Significant speedup requires reducing memory traffic (tmp elimination).',
                  transform=ax_table.transAxes, ha='center', fontsize=8, color='gray')

    plt.tight_layout()
    plt.savefig(os.path.join(BASE, 'roofline_opt_summary.png'), dpi=200, bbox_inches='tight')
    print("Saved: roofline_opt_summary.png")
    plt.close()


def main():
    print("Loading original roofline data...")
    orig_data = load_csv(ORIG_CSV)
    print(f"  Loaded {len(orig_data)} entries")

    print("Loading optimized roofline data...")
    opt_data = load_csv(OPT_CSV)
    print(f"  Loaded {len(opt_data)} entries")

    plot_roofline_comparison(orig_data, opt_data)
    plot_time_comparison(orig_data, opt_data)
    plot_bandwidth_comparison(orig_data, opt_data)
    plot_full_roofline_with_summary(orig_data, opt_data)


if __name__ == '__main__':
    main()
