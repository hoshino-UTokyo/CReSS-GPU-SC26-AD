#!/usr/bin/env python3
"""Roofline model visualization for CReSS GPU kernels on GH200."""

import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# GH200 specs
PEAK_FP32 = 59300  # GFLOPS
PEAK_FP64 = 29600  # GFLOPS
PEAK_BW = 4000     # GB/s

CSV_FILE = "roofline_summary_20260323_233343.csv"
OPT_CSV = "../../Kernel_benchmark_gpu_opt/roofline_results/roofline_summary_20260323_233343.csv"

def load_csv(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            ai = float(r['ai_dram'])
            gflops = float(r['gflops'])
            bw = float(r['bw_gbs'])
            if ai > 0 and gflops > 0:
                rows.append({
                    'benchmark': r['benchmark'],
                    'kernel': r['kernel'],
                    'ai_dram': ai,
                    'gflops': gflops,
                    'bw_gbs': bw,
                    'total_flop': int(r['total_flop']),
                    'fp64_flop': int(r['fp64_flop']),
                })
    return rows

def plot_roofline(data, title, outfile, peak_flops=PEAK_FP64):
    fig, ax = plt.subplots(figsize=(14, 9))

    # Roofline
    ai_range = np.logspace(-3, 3, 500)
    roofline = np.minimum(peak_flops, PEAK_BW * ai_range)
    ax.plot(ai_range, roofline, 'k-', linewidth=2.5, label=f'Roofline (FP64={PEAK_FP64} GFLOPS, BW={PEAK_BW} GB/s)')

    # Ridge point
    ridge_ai = peak_flops / PEAK_BW
    ax.axvline(ridge_ai, color='gray', linestyle=':', alpha=0.5)
    ax.annotate(f'Ridge point\nAI={ridge_ai:.1f}', xy=(ridge_ai, peak_flops),
                xytext=(ridge_ai * 2, peak_flops * 0.5), fontsize=8, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5))

    # Efficiency bands
    for eff, alpha in [(0.5, 0.08), (0.25, 0.06), (0.1, 0.04)]:
        roof_eff = np.minimum(peak_flops * eff, PEAK_BW * eff * ai_range)
        ax.plot(ai_range, roof_eff, '--', color='gray', alpha=0.3, linewidth=0.8)
        ax.text(ai_range[-1] * 0.7, roof_eff[-1] * 1.1, f'{int(eff*100)}%', fontsize=7, color='gray', alpha=0.5)

    # Classify kernels
    ais = [d['ai_dram'] for d in data]
    gfs = [d['gflops'] for d in data]
    sizes = [max(20, min(150, d['total_flop'] / 1e8)) for d in data]

    # Color by efficiency
    efficiencies = []
    for d in data:
        roof_val = min(peak_flops, PEAK_BW * d['ai_dram'])
        eff = d['gflops'] / roof_val if roof_val > 0 else 0
        efficiencies.append(eff)

    sc = ax.scatter(ais, gfs, c=efficiencies, cmap='RdYlGn', s=sizes,
                    alpha=0.7, edgecolors='black', linewidth=0.3,
                    vmin=0, vmax=1, zorder=5)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label('Roofline Efficiency', fontsize=10)

    # Annotate top kernels by GFLOPS
    sorted_data = sorted(zip(data, efficiencies), key=lambda x: x[0]['gflops'], reverse=True)
    annotated = []
    for d, eff in sorted_data[:15]:
        name = d['benchmark'].split('_', 1)[-1] if '_' in d['benchmark'] else d['benchmark']
        # Remove trailing suffixes like _subroutine, _s_xxx
        short = name.split('_subroutine')[0].split('_s_')[0]
        # Avoid overlapping labels
        skip = False
        for ax_ai, ax_gf in annotated:
            if abs(np.log10(d['ai_dram']) - np.log10(ax_ai)) < 0.15 and abs(np.log10(d['gflops']) - np.log10(ax_gf)) < 0.15:
                skip = True
                break
        if skip:
            continue
        ax.annotate(short, (d['ai_dram'], d['gflops']),
                    fontsize=6, alpha=0.8,
                    xytext=(5, 5), textcoords='offset points')
        annotated.append((d['ai_dram'], d['gflops']))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Arithmetic Intensity (FLOP/Byte)', fontsize=12)
    ax.set_ylabel('Performance (GFLOPS)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(1e-3, 1e2)
    ax.set_ylim(1e-1, 1e5)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(loc='upper left', fontsize=9)

    # Stats text
    n = len(data)
    mem_bound = sum(1 for e in efficiencies if e < 0.5)
    avg_eff = np.mean(efficiencies)
    stats = f'N={n} kernels | Memory-bound(<50% eff): {mem_bound} ({mem_bound/n*100:.0f}%) | Avg eff: {avg_eff:.1%}'
    ax.text(0.5, -0.08, stats, transform=ax.transAxes, fontsize=9, ha='center', color='#555')

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


def plot_bandwidth(data, outfile):
    """Bandwidth utilization histogram."""
    fig, ax = plt.subplots(figsize=(12, 6))
    bws = [d['bw_gbs'] for d in data]

    ax.hist(bws, bins=40, color='steelblue', edgecolor='white', alpha=0.8)
    ax.axvline(PEAK_BW, color='red', linestyle='--', linewidth=2, label=f'Peak BW ({PEAK_BW} GB/s)')
    ax.axvline(np.median(bws), color='orange', linestyle='--', linewidth=1.5, label=f'Median ({np.median(bws):.0f} GB/s)')

    ax.set_xlabel('DRAM Bandwidth (GB/s)', fontsize=12)
    ax.set_ylabel('Number of Kernels', fontsize=12)
    ax.set_title('CReSS GPU Kernel - DRAM Bandwidth Distribution', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


def plot_top_kernels(data, outfile, n=30):
    """Top N kernels by execution time (proxy: total_flop / gflops)."""
    for d in data:
        d['time_us'] = d['total_flop'] / (d['gflops'] * 1e3) if d['gflops'] > 0 else 0

    top = sorted(data, key=lambda x: x['time_us'], reverse=True)[:n]
    top.reverse()

    names = []
    for d in top:
        name = d['benchmark'].split('_', 1)[-1] if '_' in d['benchmark'] else d['benchmark']
        short = name.split('_subroutine')[0].split('_s_')[0]
        names.append(short)

    times = [d['time_us'] for d in top]
    effs = []
    for d in top:
        roof_val = min(PEAK_FP64, PEAK_BW * d['ai_dram'])
        effs.append(d['gflops'] / roof_val if roof_val > 0 else 0)

    colors = plt.cm.RdYlGn([e for e in effs])

    fig, ax = plt.subplots(figsize=(12, 10))
    bars = ax.barh(range(len(top)), times, color=colors, edgecolor='gray', linewidth=0.5)

    for i, (t, e) in enumerate(zip(times, effs)):
        ax.text(t + max(times) * 0.01, i, f'{e:.0%}', va='center', fontsize=7, color='#333')

    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel('Estimated Kernel Time (us)', fontsize=11)
    ax.set_title(f'Top {n} CReSS GPU Kernels by Execution Time (color = roofline efficiency)', fontsize=13)
    ax.grid(True, alpha=0.3, axis='x')

    sm = plt.cm.ScalarMappable(cmap='RdYlGn', norm=plt.Normalize(0, 1))
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label('Roofline Efficiency', fontsize=9)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


def main():
    print("Loading data...")
    data = load_csv(CSV_FILE)
    print(f"  Loaded {len(data)} kernels (with AI > 0 and GFLOPS > 0)")

    # 1. Main roofline plot
    plot_roofline(data, 'CReSS GPU Kernel Roofline Model (GH200, FP64)', 'roofline_plot.png')

    # 2. Bandwidth distribution
    plot_bandwidth(data, 'bandwidth_distribution.png')

    # 3. Top kernels by time
    plot_top_kernels(data, 'top_kernels_by_time.png')

    # 4. Summary stats
    print("\n=== Summary ===")
    ais = [d['ai_dram'] for d in data]
    gfs = [d['gflops'] for d in data]
    bws = [d['bw_gbs'] for d in data]
    print(f"  AI range: {min(ais):.4f} - {max(ais):.4f} FLOP/Byte")
    print(f"  GFLOPS range: {min(gfs):.1f} - {max(gfs):.1f}")
    print(f"  BW range: {min(bws):.0f} - {max(bws):.0f} GB/s")
    print(f"  Median AI: {np.median(ais):.4f}")
    print(f"  Median GFLOPS: {np.median(gfs):.1f}")
    print(f"  Median BW: {np.median(bws):.0f} GB/s")


if __name__ == '__main__':
    main()
