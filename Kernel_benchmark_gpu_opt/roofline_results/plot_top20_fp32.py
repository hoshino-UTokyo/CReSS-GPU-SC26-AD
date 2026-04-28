#!/usr/bin/env python3
"""
Top 10/20 kernels by measured simulation time per large timestep.
Uses actual nsys profiling data (profile_steptime, 10 large timesteps)
combined with NCU roofline data for bandwidth efficiency coloring.
"""

import csv
import sqlite3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))

# GH200 specs
PEAK_BW = 4000  # GB/s

# nsys profile ran 10 large timesteps
N_LARGE_STEPS = 10

# Initialization-only kernels to exclude
INIT_KERNELS = {'setcst3d', 'setcst4d', 'setcst2d', 'getarea', 'phycood', 'baserho'}

# Paths
NSYS_DB = os.path.join(BASE, "../../test_real/profile_steptime.sqlite")
ROOFLINE_CSV = os.path.join(BASE, "../../Kernel_benchmark_gpu/roofline_results/"
                                   "roofline_summary_20260323_233343.csv")


def extract_subroutine(kernel_name):
    """m_pgrad_s_pgrad_325_gpu -> pgrad"""
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
        SELECT s.value, COUNT(*), SUM(k.end - k.start)
        FROM CUPTI_ACTIVITY_KIND_KERNEL k
        JOIN StringIds s ON k.demangledName = s.id
        GROUP BY s.value
    """)

    sub_data = defaultdict(lambda: {'total_ns': 0, 'num_calls': 0})
    for name, calls, total_ns in cur.fetchall():
        sub = extract_subroutine(name)
        sub_data[sub]['total_ns'] += total_ns
        sub_data[sub]['num_calls'] += calls
    conn.close()

    # Exclude init kernels
    filtered = {k: v for k, v in sub_data.items() if k not in INIT_KERNELS}
    return filtered


def load_roofline_bw():
    """Load NCU roofline data to get bandwidth per subroutine."""
    bw_data = defaultdict(lambda: {'bw_gbs': [], 'ai_dram': []})
    with open(ROOFLINE_CSV) as f:
        reader = csv.DictReader(f)
        for r in reader:
            bench = r['benchmark']
            parts = bench.split('_')
            if len(parts) >= 2:
                name = '_'.join(parts[1:]).replace('_subroutine', '')
                for suffix in ['_sec1', '_sec2', '_sec3']:
                    name = name.replace(suffix, '')
                if '_s_' in name:
                    name = name.split('_s_')[1]
            else:
                name = bench

            bw = float(r['bw_gbs'])
            ai = float(r['ai_dram'])
            if bw > 0:
                bw_data[name]['bw_gbs'].append(bw)
                bw_data[name]['ai_dram'].append(ai)

    # Average per subroutine
    result = {}
    for name, data in bw_data.items():
        result[name] = {
            'bw_gbs': np.mean(data['bw_gbs']),
            'ai_dram': np.mean(data['ai_dram']),
        }
    return result


def match_roofline(sub_name, roofline):
    """Find matching roofline entry for a subroutine."""
    if sub_name in roofline:
        return roofline[sub_name]
    for rname in roofline:
        if sub_name in rname or rname in sub_name:
            return roofline[rname]
    return None


def main():
    print("Loading nsys profiling data...")
    nsys = load_nsys_data()

    print("Loading NCU roofline data...")
    roofline = load_roofline_bw()

    # Compute per-large-timestep time
    results = []
    for sub_name, data in nsys.items():
        time_per_step_ms = data['total_ns'] / N_LARGE_STEPS / 1e6  # ms
        calls_per_step = data['num_calls'] / N_LARGE_STEPS

        rf = match_roofline(sub_name, roofline)
        bw_gbs = rf['bw_gbs'] if rf else 0
        bw_eff = bw_gbs / PEAK_BW

        results.append({
            'name': sub_name,
            'time_per_step_ms': time_per_step_ms,
            'calls_per_step': calls_per_step,
            'total_ns': data['total_ns'],
            'bw_gbs': bw_gbs,
            'bw_eff': bw_eff,
        })

    results.sort(key=lambda x: x['time_per_step_ms'], reverse=True)

    total_step_ms = sum(r['time_per_step_ms'] for r in results)

    print(f"\nMeasured GPU time per large timestep: {total_step_ms:.1f} ms")
    print(f"  (from {N_LARGE_STEPS} large timesteps, init kernels excluded)")
    print(f"\nTop 20 subroutines:")
    cum = 0
    for i, r in enumerate(results[:20]):
        pct = r['time_per_step_ms'] / total_step_ms * 100
        cum += pct
        print(f"  {i+1:2d}. {r['name']:15s}  {r['time_per_step_ms']:7.2f} ms  "
              f"({pct:5.1f}%, cum {cum:5.1f}%)  "
              f"BW={r['bw_gbs']:.0f} GB/s ({r['bw_eff']:.0%})  "
              f"~{r['calls_per_step']:.0f} calls/step")

    # ================================================================
    # Plot 1: Top 10
    # ================================================================
    top10 = results[:10]
    top10_rev = top10[::-1]

    fig, ax = plt.subplots(figsize=(10, 5.5))

    names = [r['name'] for r in top10_rev]
    times = [r['time_per_step_ms'] for r in top10_rev]
    bw_effs = [r['bw_eff'] for r in top10_rev]

    colors = plt.cm.RdYlGn(bw_effs)
    bars = ax.barh(range(len(top10_rev)), times, color=colors,
                   edgecolor='gray', linewidth=0.5, height=0.7)

    ax.set_yticks(range(len(top10_rev)))
    ax.set_yticklabels(names, fontsize=11, fontweight='bold')
    ax.set_xlabel('Measured GPU time per large timestep [ms]', fontsize=12)
    ax.tick_params(axis='x', labelsize=11)
    ax.grid(True, alpha=0.3, axis='x')

    sm = plt.cm.ScalarMappable(cmap='RdYlGn', norm=plt.Normalize(0, 1))
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.02, aspect=25)
    cbar.set_label('BW efficiency', fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    plt.tight_layout()
    out1 = os.path.join(BASE, 'top10_sim_time_fp32.png')
    plt.savefig(out1, dpi=200, bbox_inches='tight')
    print(f"\nSaved: {out1}")
    plt.close()

    # ================================================================
    # Plot 2: Top 20
    # ================================================================
    top20 = results[:20]
    top20_rev = top20[::-1]

    fig, ax = plt.subplots(figsize=(12, 9))

    names = [r['name'] for r in top20_rev]
    times = [r['time_per_step_ms'] for r in top20_rev]
    bw_effs = [r['bw_eff'] for r in top20_rev]

    colors = plt.cm.RdYlGn(bw_effs)
    bars = ax.barh(range(len(top20_rev)), times, color=colors,
                   edgecolor='gray', linewidth=0.5, height=0.7)

    for i, r in enumerate(top20_rev):
        pct = r['time_per_step_ms'] / total_step_ms * 100
        ax.text(r['time_per_step_ms'] + max(times) * 0.01, i,
                f"{pct:.1f}%  BW={r['bw_gbs']:.0f} GB/s ({r['bw_eff']:.0%})",
                va='center', fontsize=7, color='#333')

    ax.set_yticks(range(len(top20_rev)))
    ax.set_yticklabels(names, fontsize=10, fontweight='bold')
    ax.set_xlabel('Measured GPU time per large timestep [ms]', fontsize=12)
    ax.set_title('CReSS GPU: Top 20 Subroutines per Large Timestep (nsys, 10 steps)',
                 fontsize=13, fontweight='bold')
    ax.tick_params(axis='x', labelsize=11)
    ax.grid(True, alpha=0.3, axis='x')

    sm = plt.cm.ScalarMappable(cmap='RdYlGn', norm=plt.Normalize(0, 1))
    cbar = plt.colorbar(sm, ax=ax, shrink=0.7, pad=0.02, aspect=25)
    cbar.set_label('BW efficiency', fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    top20_pct = sum(r['time_per_step_ms'] for r in top20) / total_step_ms * 100
    ax.text(0.98, 0.02,
            f'Top 20 = {top20_pct:.1f}% of {total_step_ms:.1f} ms/step',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    out2 = os.path.join(BASE, 'top20_sim_time_fp32.png')
    plt.savefig(out2, dpi=200, bbox_inches='tight')
    print(f"Saved: {out2}")
    plt.close()


if __name__ == '__main__':
    main()
