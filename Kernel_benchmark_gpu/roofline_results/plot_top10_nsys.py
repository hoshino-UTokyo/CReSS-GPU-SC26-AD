#!/usr/bin/env python3
"""Top-10 GPU kernel analysis based on nsys profiling data (profile_steptime)
combined with NCU roofline metrics."""

import csv
import sqlite3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from collections import defaultdict

# ── Paths ──
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
NSYS_DB = os.path.join(SCRIPT_DIR, "../../test_real/profile_steptime.sqlite")
ROOFLINE_CSV = os.path.join(SCRIPT_DIR, "roofline_summary_20260323_233343.csv")
OUT_DIR = SCRIPT_DIR

# ── GH200 specs ──
PEAK_FP64 = 29600   # GFLOPS
PEAK_BW   = 4000    # GB/s

# ── Initialization-only kernels (called only a few times, not per-step) ──
INIT_KERNELS = {'setcst3d', 'setcst4d', 'setcst2d', 'getarea', 'phycood', 'baserho'}

# ── Color palette ──
COLORS_10 = ['#2196F3', '#FF5722', '#4CAF50', '#FFC107', '#9C27B0',
             '#00BCD4', '#FF9800', '#795548', '#607D8B', '#E91E63']

def extract_subroutine(kernel_name):
    """Extract subroutine name from GPU kernel name.
    e.g., m_pgrad_s_pgrad_325_gpu -> pgrad
    """
    name = kernel_name.replace('__red', '')
    parts = name.split('_s_')
    if len(parts) >= 2:
        rest = parts[1]
        sub_parts = rest.rsplit('_', 2)
        if len(sub_parts) >= 3 and sub_parts[-1] == 'gpu':
            return sub_parts[0]
        elif len(sub_parts) >= 2:
            return sub_parts[0]
    return kernel_name


def load_nsys_data():
    """Load kernel timing from nsys sqlite, aggregate by subroutine."""
    conn = sqlite3.connect(NSYS_DB)
    cur = conn.cursor()
    cur.execute("""
        SELECT s.value, COUNT(*), SUM(k.end - k.start), AVG(k.end - k.start)
        FROM CUPTI_ACTIVITY_KIND_KERNEL k
        JOIN StringIds s ON k.demangledName = s.id
        GROUP BY s.value
        ORDER BY SUM(k.end - k.start) DESC
    """)

    sub_data = defaultdict(lambda: {'total_ns': 0, 'num_calls': 0, 'kernel_count': 0, 'kernels': []})
    for name, calls, total_ns, avg_ns in cur.fetchall():
        sub = extract_subroutine(name)
        sub_data[sub]['total_ns'] += total_ns
        sub_data[sub]['num_calls'] += calls
        sub_data[sub]['kernel_count'] += 1
        sub_data[sub]['kernels'].append({
            'name': name, 'calls': calls, 'total_ns': total_ns, 'avg_ns': avg_ns
        })
    conn.close()

    # Exclude initialization-only kernels (not called per time step)
    filtered = {k: v for k, v in sub_data.items() if k not in INIT_KERNELS}
    total_gpu = sum(v['total_ns'] for v in filtered.values())
    sorted_subs = sorted(filtered.items(), key=lambda x: x[1]['total_ns'], reverse=True)
    return sorted_subs, total_gpu


def load_roofline_data():
    """Load NCU roofline metrics, aggregate by subroutine."""
    roofline = defaultdict(lambda: {'ai_dram': [], 'gflops': [], 'bw_gbs': []})
    with open(ROOFLINE_CSV) as f:
        reader = csv.DictReader(f)
        for r in reader:
            bench = r['benchmark']
            parts = bench.split('_')
            if len(parts) >= 2:
                name = '_'.join(parts[1:]).replace('_subroutine', '')
                # Remove _secN suffix for matching
                for suffix in ['_sec1', '_sec2', '_sec3']:
                    name = name.replace(suffix, '')
                # Also remove _s_XXX pattern if present
                if '_s_' in name:
                    name = name.split('_s_')[1]
            else:
                name = bench

            ai = float(r['ai_dram'])
            gflops = float(r['gflops'])
            bw = float(r['bw_gbs'])
            if ai > 0 and gflops > 0:
                roofline[name]['ai_dram'].append(ai)
                roofline[name]['gflops'].append(gflops)
                roofline[name]['bw_gbs'].append(bw)
    return roofline


def plot_top10_time_bar(top10, total_gpu):
    """Fig 1: Top 10 kernels by total GPU time (horizontal bar chart)."""
    fig, ax = plt.subplots(figsize=(14, 7))

    names = [s[0] for s in top10][::-1]
    times_ms = [s[1]['total_ns'] / 1e6 for s in top10][::-1]
    calls = [s[1]['num_calls'] for s in top10][::-1]
    pcts = [s[1]['total_ns'] / total_gpu * 100 for s in top10][::-1]
    colors = COLORS_10[::-1]

    bars = ax.barh(range(10), times_ms, color=colors, edgecolor='white', linewidth=0.8, height=0.7)

    for i, (t, p, c) in enumerate(zip(times_ms, pcts, calls)):
        ax.text(t + max(times_ms) * 0.01, i, f'{p:.1f}%  ({c:,} calls)',
                va='center', fontsize=9, color='#333', fontweight='bold')

    ax.set_yticks(range(10))
    ax.set_yticklabels(names, fontsize=11, fontweight='bold')
    ax.set_xlabel('Total GPU Time (ms)', fontsize=12)
    ax.set_title('CReSS GPU: Top 10 Step-loop Subroutines by Total Execution Time (nsys)',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')

    # Cumulative %
    cum = 0
    top10_total = sum(s[1]['total_ns'] for s in top10) / total_gpu * 100
    ax.text(0.98, 0.02,
            f'Top 10 = {top10_total:.1f}% of step-kernel GPU time ({total_gpu/1e9:.2f} s)\n(init kernels excluded)',
            transform=ax.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    outfile = os.path.join(OUT_DIR, 'top10_time_bar.png')
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


def plot_top10_pie(top10, total_gpu):
    """Fig 2: Top 10 GPU time breakdown (pie chart)."""
    fig, ax = plt.subplots(figsize=(10, 10))

    names = [s[0] for s in top10]
    times = [s[1]['total_ns'] for s in top10]
    other = total_gpu - sum(times)
    names.append('Others')
    times.append(other)
    colors = COLORS_10 + ['#E0E0E0']

    def fmt_pct(pct):
        return f'{pct:.1f}%' if pct >= 2 else ''

    wedges, texts, autotexts = ax.pie(
        times, labels=names, colors=colors, autopct=fmt_pct,
        startangle=90, pctdistance=0.82,
        wedgeprops=dict(edgecolor='white', linewidth=1.5))

    for t in texts:
        t.set_fontsize(10)
        t.set_fontweight('bold')
    for t in autotexts:
        t.set_fontsize(8)

    ax.set_title('CReSS GPU: Step-loop Kernel Time Distribution (nsys, init excluded)',
                 fontsize=14, fontweight='bold', pad=20)

    # Add legend with times
    legend_labels = []
    for i, (n, t) in enumerate(zip(names, times)):
        pct = t / total_gpu * 100
        legend_labels.append(f'{n}: {t/1e6:.0f} ms ({pct:.1f}%)')
    ax.legend(wedges, legend_labels, loc='center left', bbox_to_anchor=(1.0, 0.5),
              fontsize=9, framealpha=0.9)

    plt.tight_layout()
    outfile = os.path.join(OUT_DIR, 'top10_time_pie.png')
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


def plot_top10_roofline(top10, roofline):
    """Fig 3: Roofline plot highlighting top-10 subroutines."""
    fig, ax = plt.subplots(figsize=(14, 9))

    # Draw roofline
    ai_range = np.logspace(-3, 3, 500)
    roofline_curve = np.minimum(PEAK_FP64, PEAK_BW * ai_range)
    ax.plot(ai_range, roofline_curve, 'k-', linewidth=2.5,
            label=f'Roofline (FP64={PEAK_FP64} GFLOPS, BW={PEAK_BW} GB/s)')

    # Ridge point
    ridge_ai = PEAK_FP64 / PEAK_BW
    ax.axvline(ridge_ai, color='gray', linestyle=':', alpha=0.5)
    ax.annotate(f'Ridge AI={ridge_ai:.1f}', xy=(ridge_ai, PEAK_FP64),
                xytext=(ridge_ai * 3, PEAK_FP64 * 0.4), fontsize=9, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5))

    # Efficiency lines
    for eff in [0.5, 0.25, 0.1]:
        roof_eff = np.minimum(PEAK_FP64 * eff, PEAK_BW * eff * ai_range)
        ax.plot(ai_range, roof_eff, '--', color='gray', alpha=0.25, linewidth=0.8)
        ax.text(ai_range[-1] * 0.7, roof_eff[-1] * 1.1, f'{int(eff*100)}%',
                fontsize=7, color='gray', alpha=0.5)

    # Plot top-10 subroutines
    top10_names = [s[0] for s in top10]
    for i, (sub_name, sub_info) in enumerate(top10):
        # Try to find matching roofline data
        matched = None
        for rname in roofline:
            if sub_name == rname or sub_name in rname or rname in sub_name:
                matched = rname
                break

        if matched and roofline[matched]['ai_dram']:
            ais = roofline[matched]['ai_dram']
            gfs = roofline[matched]['gflops']
            # Plot individual kernel data points
            ax.scatter(ais, gfs, color=COLORS_10[i], s=100, alpha=0.8,
                       edgecolors='black', linewidth=0.8, zorder=10,
                       label=f'{sub_name} ({sub_info["total_ns"]/1e6:.0f} ms)')
            # Annotate with name
            avg_ai = np.mean(ais)
            avg_gf = np.mean(gfs)
            ax.annotate(sub_name, (avg_ai, avg_gf), fontsize=8, fontweight='bold',
                        color=COLORS_10[i],
                        xytext=(8, 8), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Arithmetic Intensity (FLOP/Byte)', fontsize=12)
    ax.set_ylabel('Performance (GFLOPS)', fontsize=12)
    ax.set_title('CReSS GPU: Top 10 Step-loop Subroutines on Roofline Model (GH200, FP64)',
                 fontsize=14, fontweight='bold')
    ax.set_xlim(1e-2, 1e2)
    ax.set_ylim(1e1, 5e4)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(loc='upper left', fontsize=8, ncol=1, framealpha=0.9)

    plt.tight_layout()
    outfile = os.path.join(OUT_DIR, 'top10_roofline.png')
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


def plot_top10_detail_table(top10, roofline, total_gpu):
    """Fig 4: Summary table + per-call time comparison."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={'width_ratios': [1.2, 1]})

    # Left: stacked bar showing avg time per call vs num calls
    names = [s[0] for s in top10]
    avg_us = [s[1]['total_ns'] / s[1]['num_calls'] / 1e3 for s in top10]  # microseconds
    num_calls = [s[1]['num_calls'] for s in top10]

    # Bar chart: avg time per call
    bars = ax1.barh(range(10), avg_us, color=COLORS_10, edgecolor='white', linewidth=0.8, height=0.7)
    for i, (t, n) in enumerate(zip(avg_us, num_calls)):
        ax1.text(t + max(avg_us) * 0.01, i,
                 f'{t:.0f} us x {n:,} calls',
                 va='center', fontsize=8, color='#333')

    ax1.set_yticks(range(10))
    ax1.set_yticklabels(names, fontsize=10, fontweight='bold')
    ax1.set_xlabel('Average Kernel Time per Call (us)', fontsize=11)
    ax1.set_title('Avg. Time per Call', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='x')

    # Right: Roofline efficiency + bandwidth
    efficiencies = []
    bandwidths = []
    for sub_name, sub_info in top10:
        matched = None
        for rname in roofline:
            if sub_name == rname or sub_name in rname or rname in sub_name:
                matched = rname
                break
        if matched and roofline[matched]['ai_dram']:
            avg_ai = np.mean(roofline[matched]['ai_dram'])
            avg_gf = np.mean(roofline[matched]['gflops'])
            avg_bw = np.mean(roofline[matched]['bw_gbs'])
            roof_val = min(PEAK_FP64, PEAK_BW * avg_ai)
            eff = avg_gf / roof_val if roof_val > 0 else 0
            efficiencies.append(eff)
            bandwidths.append(avg_bw)
        else:
            efficiencies.append(0)
            bandwidths.append(0)

    x = np.arange(10)
    width = 0.35
    bars1 = ax2.bar(x - width/2, [e * 100 for e in efficiencies], width,
                     color=COLORS_10, edgecolor='white', linewidth=0.8, label='Roofline Eff. (%)')

    ax2_twin = ax2.twinx()
    bars2 = ax2_twin.bar(x + width/2, bandwidths, width,
                          color=[c + '80' for c in COLORS_10] if False else
                          [plt.matplotlib.colors.to_rgba(c, 0.4) for c in COLORS_10],
                          edgecolor='gray', linewidth=0.5, label='DRAM BW (GB/s)')

    ax2.set_xticks(x)
    ax2.set_xticklabels(names, fontsize=9, rotation=45, ha='right', fontweight='bold')
    ax2.set_ylabel('Roofline Efficiency (%)', fontsize=11)
    ax2_twin.set_ylabel('DRAM Bandwidth (GB/s)', fontsize=11)
    ax2.set_title('Efficiency & Bandwidth', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 100)
    ax2_twin.set_ylim(0, PEAK_BW)
    ax2.axhline(50, color='orange', linestyle='--', alpha=0.5, linewidth=0.8)
    ax2_twin.axhline(PEAK_BW, color='red', linestyle='--', alpha=0.3, linewidth=0.8)
    ax2.grid(True, alpha=0.2, axis='y')

    # Legends
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=8)

    fig.suptitle('CReSS GPU: Top 10 Kernel Detailed Analysis (nsys + NCU)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    outfile = os.path.join(OUT_DIR, 'top10_detail.png')
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


def plot_top10_cumulative(sorted_subs, total_gpu):
    """Fig 5: Cumulative GPU time coverage (Pareto chart)."""
    fig, ax1 = plt.subplots(figsize=(14, 7))

    n_show = min(30, len(sorted_subs))
    names = [s[0] for s in sorted_subs[:n_show]]
    times_ms = [s[1]['total_ns'] / 1e6 for s in sorted_subs[:n_show]]
    pcts = [s[1]['total_ns'] / total_gpu * 100 for s in sorted_subs[:n_show]]

    # Bar chart
    colors_bar = COLORS_10 + ['#B0BEC5'] * (n_show - 10)
    bars = ax1.bar(range(n_show), times_ms, color=colors_bar, edgecolor='white', linewidth=0.5)

    # Cumulative line
    ax2 = ax1.twinx()
    cum_pct = np.cumsum(pcts)
    ax2.plot(range(n_show), cum_pct, 'r-o', markersize=4, linewidth=2, label='Cumulative %')
    ax2.axhline(80, color='red', linestyle='--', alpha=0.4, linewidth=1)
    ax2.text(n_show - 1, 81, '80%', fontsize=9, color='red', alpha=0.6)
    ax2.axhline(90, color='darkred', linestyle='--', alpha=0.3, linewidth=1)
    ax2.text(n_show - 1, 91, '90%', fontsize=9, color='darkred', alpha=0.5)

    # Mark top 10 boundary
    ax1.axvline(9.5, color='blue', linestyle=':', alpha=0.5, linewidth=1.5)
    ax1.text(9.5, max(times_ms) * 0.95, f'  Top 10: {cum_pct[9]:.1f}%',
             fontsize=9, color='blue', fontweight='bold')

    ax1.set_xticks(range(n_show))
    ax1.set_xticklabels(names, fontsize=8, rotation=60, ha='right')
    ax1.set_ylabel('Total GPU Time (ms)', fontsize=12)
    ax2.set_ylabel('Cumulative GPU Time (%)', fontsize=12, color='red')
    ax2.set_ylim(0, 105)
    ax1.set_title('CReSS GPU: Step-loop Kernel Time Pareto Chart (Top 30, init excluded)',
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.2, axis='y')
    ax2.legend(loc='center right', fontsize=9)

    plt.tight_layout()
    outfile = os.path.join(OUT_DIR, 'top10_pareto.png')
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


def main():
    print("Loading nsys profiling data...")
    sorted_subs, total_gpu = load_nsys_data()
    print(f"  Total GPU time: {total_gpu/1e9:.3f} s")
    print(f"  Unique subroutines: {len(sorted_subs)}")

    top10 = sorted_subs[:10]
    print("\nTop 10 subroutines:")
    cum = 0
    for i, (name, data) in enumerate(top10):
        pct = data['total_ns'] / total_gpu * 100
        cum += pct
        print(f"  {i+1:2d}. {name:<15} {data['total_ns']/1e6:8.1f} ms  "
              f"({data['num_calls']:6d} calls, {data['kernel_count']:2d} kernels)  "
              f"{pct:.1f}% (cum {cum:.1f}%)")

    print("\nLoading roofline data...")
    roofline = load_roofline_data()
    print(f"  Loaded {len(roofline)} subroutines with roofline metrics")

    print("\nGenerating plots...")
    plot_top10_time_bar(top10, total_gpu)
    plot_top10_pie(top10, total_gpu)
    plot_top10_roofline(top10, roofline)
    plot_top10_detail_table(top10, roofline, total_gpu)
    plot_top10_cumulative(sorted_subs, total_gpu)

    print("\nDone! Generated 5 plots:")
    print("  1. top10_time_bar.png   - Top 10 by total GPU time")
    print("  2. top10_time_pie.png   - GPU time distribution pie chart")
    print("  3. top10_roofline.png   - Top 10 on roofline model")
    print("  4. top10_detail.png     - Detailed: avg time, efficiency, bandwidth")
    print("  5. top10_pareto.png     - Pareto chart (top 30)")


if __name__ == '__main__':
    main()
