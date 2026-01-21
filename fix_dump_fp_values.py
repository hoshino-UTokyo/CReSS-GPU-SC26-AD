#!/usr/bin/env python3
"""
Fix dump code to dump actual values instead of fp* index values.

Pattern: call dump_scalar_i('fpXXX', fpXXX) -> call dump_scalar_i('XXX', XXX) or dump_scalar_r('XXX', XXX)

The fp* variables are indices for getiname/getrname functions.
After getiname(fpXXX, XXX), the actual value is in XXX.

Real values (use dump_scalar_r):
- dxiv, dyiv, dziv, dx, dy, dz, and other grid spacing related
- Any variable obtained via getrname

Integer values (use dump_scalar_i):
- Options like advopt, mpopt, mfcopt, trnopt, etc.
- Boundary conditions like wbc, ebc, sbc, nbc
- Any variable obtained via getiname
"""

import re
import os
import sys

# Real variables (obtained via getrname)
REAL_VARS = {
    'dxiv', 'dyiv', 'dziv', 'dx', 'dy', 'dz',
    'dtb', 'dts', 'dtl', 'dtr',
    'thresq', 'dtcmph', 'dtvcul', 'dtgrd', 'dtsfc',
    'zsfc', 'zflat', 'stime', 'etime',
    'gpvitv', 'gsmitv', 'aslitv', 'rdritv', 'sstint',
    'dmpitv', 'monitv', 'resitv',
    'alpha1', 'alpha2', 'vspgpv', 'nggdmp', 'lspvar', 'vspbar', 'botgpv',
    # Add more as needed
}

def get_actual_var_name(fp_var):
    """Extract actual variable name from fp* variable name."""
    if fp_var.startswith('fp'):
        return fp_var[2:]
    return fp_var

def is_real_var(var_name):
    """Check if variable is a real type (should use dump_scalar_r)."""
    return var_name.lower() in REAL_VARS

def fix_dump_line(line):
    """Fix a single dump_scalar line if it dumps fp* value."""
    # Pattern: call dump_scalar_i('fpXXX', fpXXX)
    # or: call dump_scalar_i('fpXXX',fpXXX)
    pattern = r"call\s+dump_scalar_i\s*\(\s*'(fp[a-zA-Z0-9_]+)'\s*,\s*(fp[a-zA-Z0-9_]+)\s*\)"

    match = re.search(pattern, line, re.IGNORECASE)
    if match:
        fp_name = match.group(1)
        fp_var = match.group(2)

        # Get actual variable name
        actual_name = get_actual_var_name(fp_name)
        actual_var = get_actual_var_name(fp_var)

        # Determine if it should be dump_scalar_r or dump_scalar_i
        if is_real_var(actual_name):
            new_call = f"call dump_scalar_r('{actual_name}', {actual_var})"
        else:
            new_call = f"call dump_scalar_i('{actual_name}', {actual_var})"

        # Replace in line
        new_line = re.sub(pattern, new_call, line, flags=re.IGNORECASE)
        return new_line, True

    return line, False

def fix_file(filepath):
    """Fix all fp* dump calls in a file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    modified = False
    new_lines = []
    changes = []

    for i, line in enumerate(lines):
        new_line, changed = fix_dump_line(line)
        new_lines.append(new_line)
        if changed:
            modified = True
            changes.append((i+1, line.strip(), new_line.strip()))

    if modified:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        return changes

    return []

def main():
    src_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Src')

    if not os.path.isdir(src_dir):
        print(f"Error: {src_dir} not found")
        sys.exit(1)

    total_changes = 0
    files_modified = 0

    for filename in sorted(os.listdir(src_dir)):
        if filename.endswith('.f90'):
            filepath = os.path.join(src_dir, filename)
            changes = fix_file(filepath)
            if changes:
                files_modified += 1
                print(f"\n{filename}:")
                for line_num, old, new in changes:
                    print(f"  Line {line_num}:")
                    print(f"    - {old}")
                    print(f"    + {new}")
                    total_changes += 1

    print(f"\n{'='*60}")
    print(f"Total: {total_changes} changes in {files_modified} files")

if __name__ == '__main__':
    main()
