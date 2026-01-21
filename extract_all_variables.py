#!/usr/bin/env python3
"""
Extract variable lists for all kernels using default(none) method.
This script:
1. Finds all source files with OpenMP parallel sections
2. Temporarily changes default(shared) to default(none)
3. Compiles to find missing variable errors
4. Iteratively adds variables until compilation succeeds
5. Creates variable_list.txt in the corresponding Kernel_benchmark directory
6. Restores the original source file
"""

import os
import re
import subprocess
import shutil
from pathlib import Path

SRC_DIR = Path("/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload/Src")
KERNEL_DIR = Path("/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload/Kernel_benchmark")
COMPILER = "nvfortran"
COMPILER_FLAGS = ["-fast", "-mp", "-Mbyteswapio", "-mcmodel=medium", "-c"]

def find_omp_parallel_sections(source_file):
    """Find all !$omp parallel lines in the source file."""
    with open(source_file, 'r') as f:
        content = f.read()

    # Find all !$omp parallel default(shared) patterns
    pattern = r'!\$omp\s+parallel\s+default\(shared\)'
    matches = list(re.finditer(pattern, content, re.IGNORECASE))
    return matches

def get_missing_variables(source_file):
    """Compile and extract missing variable names from errors."""
    result = subprocess.run(
        [COMPILER] + COMPILER_FLAGS + [str(source_file)],
        capture_output=True,
        text=True,
        cwd=SRC_DIR
    )

    # Parse NVFORTRAN-S-0155 errors
    errors = result.stderr + result.stdout
    pattern = r'NVFORTRAN-S-0155-(\w+)\s+must appear'
    matches = re.findall(pattern, errors)
    return list(set(matches))

def update_omp_directive(source_file, shared_vars):
    """Update the !$omp parallel directive with shared variables."""
    with open(source_file, 'r') as f:
        content = f.read()

    # Build the shared clause
    if shared_vars:
        shared_clause = "shared(" + ",".join(shared_vars) + ")"
    else:
        shared_clause = ""

    # Replace pattern - handle both single line and continued lines
    # First, try to match existing default(none) with shared clause
    pattern = r'!\$omp\s+parallel\s+default\(none\)(\s*&\s*\n!\$omp&\s*shared\([^)]+\))?'

    if shared_vars:
        replacement = f"!$omp parallel default(none) {shared_clause}"
    else:
        replacement = "!$omp parallel default(none)"

    new_content = re.sub(pattern, replacement, content, flags=re.IGNORECASE)

    # Also handle the case where it's still default(shared)
    pattern2 = r'!\$omp\s+parallel\s+default\(shared\)(\s+private\([^)]+\))?'
    def replace_func(m):
        private_clause = m.group(1) if m.group(1) else ""
        if shared_vars:
            return f"!$omp parallel default(none){private_clause} {shared_clause}"
        else:
            return f"!$omp parallel default(none){private_clause}"

    new_content = re.sub(pattern2, replace_func, new_content, flags=re.IGNORECASE)

    with open(source_file, 'w') as f:
        f.write(new_content)

def process_kernel(source_file, kernel_dir):
    """Process a single kernel file."""
    print(f"\nProcessing: {source_file.name}")

    # Check if file has OpenMP parallel sections
    matches = find_omp_parallel_sections(source_file)
    if not matches:
        print(f"  No !$omp parallel default(shared) found, skipping")
        return None

    # Backup original file
    backup_file = source_file.with_suffix('.f90.bak')
    shutil.copy(source_file, backup_file)

    try:
        # Change to default(none)
        with open(source_file, 'r') as f:
            content = f.read()

        # Replace default(shared) with default(none), preserving private clause
        new_content = re.sub(
            r'(!\$omp\s+parallel\s+)default\(shared\)',
            r'\1default(none)',
            content,
            flags=re.IGNORECASE
        )

        with open(source_file, 'w') as f:
            f.write(new_content)

        # Iteratively find and add variables
        all_vars = set()
        max_iterations = 50

        for iteration in range(max_iterations):
            missing = get_missing_variables(source_file)
            if not missing:
                break

            all_vars.update(missing)
            update_omp_directive(source_file, all_vars)
            print(f"  Iteration {iteration+1}: found {len(missing)} vars, total: {len(all_vars)}")

        # Create variable list file
        var_list_file = kernel_dir / "variable_list.txt"
        with open(var_list_file, 'w') as f:
            f.write(f"# Variable List for {source_file.name}\n")
            f.write(f"# Discovered via default(none) method\n")
            f.write(f"# Source: Src/{source_file.name}\n\n")
            f.write("## Shared variables (to be dumped)\n")
            for var in sorted(all_vars):
                f.write(f"{var}\n")

        print(f"  Created: {var_list_file}")
        print(f"  Total variables: {len(all_vars)}")

        return all_vars

    finally:
        # Restore original file
        shutil.move(backup_file, source_file)
        print(f"  Restored original file")

def find_kernel_dir(source_file):
    """Find the corresponding Kernel_benchmark directory for a source file."""
    base_name = source_file.stem  # e.g., "baserho"

    # Search for matching directory
    for d in KERNEL_DIR.iterdir():
        if d.is_dir() and base_name in d.name.lower():
            return d

    return None

def main():
    # Get list of source files to process
    source_files = sorted(SRC_DIR.glob("*.f90"))

    processed = 0
    skipped = 0

    for source_file in source_files:
        # Find corresponding kernel directory
        kernel_dir = find_kernel_dir(source_file)

        if kernel_dir is None:
            continue

        # Check if variable_list.txt already exists
        var_list_file = kernel_dir / "variable_list.txt"
        if var_list_file.exists():
            print(f"Skipping {source_file.name}: variable_list.txt already exists")
            skipped += 1
            continue

        result = process_kernel(source_file, kernel_dir)
        if result is not None:
            processed += 1

    print(f"\n=== Summary ===")
    print(f"Processed: {processed}")
    print(f"Skipped: {skipped}")

if __name__ == "__main__":
    main()
