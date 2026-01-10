#!/usr/bin/env python3
"""
Check dump coverage by comparing variable_list.txt with actual dump calls in source files.
"""

import os
import re
from pathlib import Path

SRC_DIR = Path("/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload/Src")
KERNEL_DIR = Path("/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload/Kernel_benchmark")

def get_required_vars(kernel_dir):
    """Read variable_list.txt and extract required variables."""
    var_file = kernel_dir / "variable_list.txt"
    if not var_file.exists():
        return set()

    vars = set()
    with open(var_file, 'r') as f:
        in_vars_section = False
        for line in f:
            line = line.strip()
            if line.startswith("## Shared variables"):
                in_vars_section = True
                continue
            if line.startswith("##"):
                in_vars_section = False
                continue
            if in_vars_section and line and not line.startswith("#"):
                # Extract variable name (first word)
                var = line.split()[0] if line.split() else ""
                if var and var.isidentifier():
                    vars.add(var.lower())
    return vars

def get_dumped_vars(source_file):
    """Extract variables that are being dumped from source file."""
    if not source_file.exists():
        return set()

    with open(source_file, 'r') as f:
        content = f.read()

    dumped = set()

    # Find dump_scalar_i calls
    for match in re.finditer(r"call\s+dump_scalar_i\s*\(\s*'([^']+)'", content, re.IGNORECASE):
        dumped.add(match.group(1).lower())

    # Find dump_scalar_r calls
    for match in re.finditer(r"call\s+dump_scalar_r\s*\(\s*'([^']+)'", content, re.IGNORECASE):
        dumped.add(match.group(1).lower())

    # Find dump_scalar_s calls (string)
    for match in re.finditer(r"call\s+dump_scalar_s\s*\(\s*'([^']+)'", content, re.IGNORECASE):
        dumped.add(match.group(1).lower())

    # Find dump_scalar_c calls (character)
    for match in re.finditer(r"call\s+dump_scalar_c\s*\(\s*'([^']+)'", content, re.IGNORECASE):
        dumped.add(match.group(1).lower())

    # Find dump_array_* calls - extract the filename and map to variable
    # Pattern matches dump_array_2d, dump_array_3d, dump_array_2d_i, dump_array_3d_i, dump_array_2d_int, etc.
    for match in re.finditer(r"call\s+dump_array_\d+d(?:_[ir]|_int)?\s*\(\s*'([^']+)'", content, re.IGNORECASE):
        filename = match.group(1).lower()
        # Remove .bin extension and _in/_ref suffix
        var = filename.replace('.bin', '').replace('_in', '').replace('_ref', '')
        dumped.add(var)

    return dumped

def find_source_file(kernel_name):
    """Find source file for a kernel."""
    # Extract base name from kernel directory
    # e.g., "044_buoywb_s_buoywb" -> "buoywb"
    parts = kernel_name.split('_')

    # Try different patterns
    for i, part in enumerate(parts):
        if part.startswith('s') and i > 0:
            # Found "s_xxx" pattern, use the name before
            base = parts[i-1]
            source = SRC_DIR / f"{base}.f90"
            if source.exists():
                return source

    # Try subroutine pattern
    for part in parts:
        if part == "subroutine":
            idx = parts.index("subroutine")
            if idx > 1:
                base = parts[idx-1]
                source = SRC_DIR / f"{base}.f90"
                if source.exists():
                    return source

    # Fallback: try each part
    for part in parts[1:]:  # Skip number prefix
        if part not in ['s', 'subroutine', 'sec1', 'sec2', 'sec3', 'sec4', 'sec5', 'sec6', 'sec7', 'sec8', 'r8']:
            source = SRC_DIR / f"{part}.f90"
            if source.exists():
                return source

    return None

def main():
    missing_report = []

    for kernel_dir in sorted(KERNEL_DIR.iterdir()):
        if not kernel_dir.is_dir():
            continue

        kernel_name = kernel_dir.name
        required = get_required_vars(kernel_dir)

        if not required:
            continue

        source_file = find_source_file(kernel_name)
        if not source_file:
            continue

        dumped = get_dumped_vars(source_file)

        # Find missing variables
        missing = required - dumped

        # Filter out known non-dumpable items
        skip_vars = {'i', 'j', 'k'}  # loop indices
        missing = missing - skip_vars

        if missing:
            missing_report.append({
                'kernel': kernel_name,
                'source': source_file.name,
                'required': len(required),
                'dumped': len(dumped),
                'missing': sorted(missing)
            })

    # Print report
    print("=" * 80)
    print("DUMP COVERAGE REPORT")
    print("=" * 80)
    print(f"\nTotal kernels checked: {len(list(KERNEL_DIR.iterdir()))}")
    print(f"Kernels with missing dumps: {len(missing_report)}")
    print()

    for item in missing_report:
        print(f"\n{item['kernel']} ({item['source']})")
        print(f"  Required: {item['required']}, Dumped: {item['dumped']}")
        print(f"  Missing: {', '.join(item['missing'])}")

    return missing_report

if __name__ == "__main__":
    main()
