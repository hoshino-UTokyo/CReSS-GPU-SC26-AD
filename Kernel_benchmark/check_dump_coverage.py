#!/usr/bin/env python3
"""
Check if all variables in variable_list.txt are dumped in the source file.
"""

import os
import re
import glob

def parse_variable_list(filepath):
    """Parse variable_list.txt and extract shared variables."""
    variables = []
    in_shared_section = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if '## Shared variables' in line:
                in_shared_section = True
                continue
            if line.startswith('##') and in_shared_section:
                break  # End of shared section
            if in_shared_section and line and not line.startswith('#'):
                # Extract variable name (first word)
                var = line.split()[0] if line.split() else ''
                if var and not var.startswith('Total'):
                    variables.append(var.lower())
    return variables

def has_dump_section(content):
    """Check if source file has dump instrumentation."""
    return 'dump_call_count' in content.lower() or 'call dump_init' in content.lower()

def parse_source_dump(filepath, kernel_name):
    """Parse source file and extract dumped variables."""
    dumped_vars = set()

    with open(filepath, 'r') as f:
        content = f.read()

    # Check if dump section exists
    has_dump = has_dump_section(content)

    # Pattern: dump_scalar_i('varname', varname) or dump_scalar_r/s/c('varname', varname)
    pattern1 = r"call\s+dump_scalar_[ircs]\s*\(\s*'(\w+)'\s*,\s*(\w+)"
    for m in re.finditer(pattern1, content, re.IGNORECASE):
        # Use the actual variable name (second group), not the string label
        dumped_vars.add(m.group(2).lower())

    # Pattern: dump_array_*d('varname.bin', varname, ...) - also matches _i, _int suffixes
    pattern2 = r"call\s+dump_array_\d+d(?:_i(?:nt)?)?\s*\(\s*'[^']+'\s*,\s*(\w+)"
    for m in re.finditer(pattern2, content, re.IGNORECASE):
        dumped_vars.add(m.group(1).lower())

    return list(dumped_vars), has_dump

def find_source_file(src_dir, kernel_dir_name):
    """Find the source file for a kernel."""
    # Parse kernel directory name: XXX_name_subroutine or XXX_name_s_name
    parts = kernel_dir_name.split('_', 1)
    if len(parts) < 2:
        return None

    name_part = parts[1]

    # Try common patterns
    possible_names = []

    # Extract base name
    if '_subroutine' in name_part:
        base = name_part.replace('_subroutine', '').replace('_sec1', '').replace('_sec2', '')
        possible_names.append(base + '.f90')
    elif '_s_' in name_part:
        base = name_part.split('_s_')[0]
        possible_names.append(base + '.f90')
    else:
        base = name_part.replace('_sec1', '').replace('_sec2', '')
        possible_names.append(base + '.f90')

    for name in possible_names:
        path = os.path.join(src_dir, name)
        if os.path.exists(path):
            return path

    return None

def main():
    base_dir = '/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload'
    kernel_dir = os.path.join(base_dir, 'Kernel_benchmark')
    src_dir = os.path.join(base_dir, 'Src')

    # Find all variable_list.txt files
    var_lists = glob.glob(os.path.join(kernel_dir, '*/variable_list.txt'))
    var_lists.sort()

    print(f"Found {len(var_lists)} kernel directories with variable_list.txt\n")

    no_dump = []       # Files with no dump section
    incomplete = []    # Files with incomplete dump
    complete = []      # Files with complete dump
    not_found = []     # Source files not found

    for var_list_path in var_lists:
        kernel_dir_name = os.path.basename(os.path.dirname(var_list_path))

        # Parse variable list
        required_vars = parse_variable_list(var_list_path)
        if not required_vars:
            continue

        # Find source file
        src_file = find_source_file(src_dir, kernel_dir_name)
        if not src_file:
            not_found.append((kernel_dir_name, required_vars))
            continue

        # Parse dumped variables
        dumped_vars, has_dump = parse_source_dump(src_file, kernel_dir_name)

        if not has_dump:
            no_dump.append((kernel_dir_name, os.path.basename(src_file), required_vars))
            continue

        # Find missing variables
        missing = [v for v in required_vars if v not in dumped_vars]

        if missing:
            incomplete.append((kernel_dir_name, os.path.basename(src_file), missing, required_vars, dumped_vars))
        else:
            complete.append((kernel_dir_name, os.path.basename(src_file)))

    # Output results
    print("=" * 70)
    print(f"COMPLETE: {len(complete)} kernels have full dump coverage")
    print("=" * 70)

    if no_dump:
        print(f"\n{'='*70}")
        print(f"NO DUMP SECTION: {len(no_dump)} kernels (need dump instrumentation)")
        print("=" * 70)
        for kernel, src, required in no_dump[:20]:  # Show first 20
            print(f"  {kernel} ({src}): {len(required)} vars needed")
        if len(no_dump) > 20:
            print(f"  ... and {len(no_dump) - 20} more")

    if incomplete:
        print(f"\n{'='*70}")
        print(f"INCOMPLETE DUMP: {len(incomplete)} kernels (missing some variables)")
        print("=" * 70)
        total_missing = 0
        for kernel, src, missing, required, dumped in incomplete:
            print(f"\n{kernel} ({src}):")
            print(f"  Required: {len(required)}, Dumped: {len(dumped)}, Missing: {len(missing)}")
            print(f"  Missing: {', '.join(sorted(missing))}")
            total_missing += len(missing)
        print(f"\nTotal missing variables: {total_missing}")

    if not_found:
        print(f"\n{'='*70}")
        print(f"SOURCE NOT FOUND: {len(not_found)} kernels")
        print("=" * 70)
        for kernel, required in not_found[:10]:
            print(f"  {kernel}")
        if len(not_found) > 10:
            print(f"  ... and {len(not_found) - 10} more")

    print(f"\n{'='*70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"  Complete:     {len(complete)}")
    print(f"  Incomplete:   {len(incomplete)}")
    print(f"  No dump:      {len(no_dump)}")
    print(f"  Not found:    {len(not_found)}")
    print(f"  Total:        {len(var_lists)}")

    return len(incomplete) + len(no_dump)

if __name__ == '__main__':
    exit(main())
