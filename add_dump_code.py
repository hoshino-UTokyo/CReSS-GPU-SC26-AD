#!/usr/bin/env python3
"""
Script to add dump code to OpenMP parallel sections for benchmark data generation.
Reads call counts from omp_profile.txt and adds dump instrumentation to source files.
"""

import os
import re
import sys
from pathlib import Path
from collections import defaultdict

def parse_omp_profile(profile_path):
    """Parse omp_profile.txt and return dict of filename -> (subroutine, count)."""
    profiles = {}
    with open(profile_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        # Skip header lines
        if line.startswith('=') or line.startswith('-') or 'ID' in line or not line.strip():
            continue

        parts = line.split()
        if len(parts) >= 4:
            try:
                profile_id = int(parts[0])
                filename = parts[1]
                subroutine = parts[2]
                count = int(parts[3])

                if filename not in profiles:
                    profiles[filename] = []
                profiles[filename].append({
                    'id': profile_id,
                    'subroutine': subroutine,
                    'count': count
                })
            except (ValueError, IndexError):
                continue

    return profiles

def parse_subroutine_args(content, subroutine_name):
    """Parse subroutine arguments and their intent/dimensions."""
    args = {
        'scalars_int': [],
        'scalars_real': [],
        'arrays_2d': [],
        'arrays_3d': [],
        'arrays_4d': [],
        'input_arrays': [],
        'output_arrays': [],
        'inout_arrays': []
    }

    # Find subroutine declaration
    subr_pattern = rf'subroutine\s+{subroutine_name}\s*\((.*?)\)'
    match = re.search(subr_pattern, content, re.IGNORECASE | re.DOTALL)
    if not match:
        return args

    # Get argument names
    arg_list = match.group(1).replace('&', '').replace('\n', ' ')
    arg_names = [a.strip() for a in arg_list.split(',') if a.strip()]

    # Find declarations for each argument
    for arg in arg_names:
        # Look for intent declarations
        # Pattern: real, intent(in) :: varname(...)
        intent_pattern = rf'(integer|real).*intent\s*\(\s*(in|out|inout)\s*\).*::\s*{arg}\s*(\([^)]+\))?'
        match = re.search(intent_pattern, content, re.IGNORECASE)

        if match:
            var_type = match.group(1).lower()
            intent = match.group(2).lower()
            dims = match.group(3)

            # Determine dimensions
            if dims:
                dim_count = dims.count(',') + 1
            else:
                dim_count = 0

            info = {
                'name': arg,
                'type': var_type,
                'intent': intent,
                'dims': dim_count,
                'dim_str': dims if dims else ''
            }

            if dim_count == 0:
                if var_type == 'integer':
                    args['scalars_int'].append(info)
                else:
                    args['scalars_real'].append(info)
            elif dim_count == 2:
                args['arrays_2d'].append(info)
            elif dim_count == 3:
                args['arrays_3d'].append(info)
            elif dim_count == 4:
                args['arrays_4d'].append(info)

            if intent == 'in':
                args['input_arrays'].append(info)
            elif intent == 'out':
                args['output_arrays'].append(info)
            elif intent == 'inout':
                args['inout_arrays'].append(info)

    return args

def extract_array_bounds(dim_str):
    """Extract array bounds from dimension string like (0:ni+1,0:nj+1,1:nk)."""
    if not dim_str:
        return []

    # Remove parentheses
    dim_str = dim_str.strip('()')
    dims = dim_str.split(',')

    bounds = []
    for dim in dims:
        dim = dim.strip()
        if ':' in dim:
            parts = dim.split(':')
            bounds.append((parts[0].strip(), parts[1].strip()))
        else:
            bounds.append(('1', dim.strip()))

    return bounds

def generate_dump_calls(args, kernel_name):
    """Generate dump call statements for the given arguments."""
    input_dumps = []
    output_dumps = []

    # Dump scalars
    for s in args['scalars_int']:
        input_dumps.append(f"  call dump_scalar_i('{s['name']}', {s['name']})")
    for s in args['scalars_real']:
        input_dumps.append(f"  call dump_scalar_r('{s['name']}', {s['name']})")

    # Dump input arrays
    for arr in args['input_arrays']:
        if arr['dims'] > 0:
            bounds = extract_array_bounds(arr['dim_str'])
            if arr['dims'] == 2 and len(bounds) == 2:
                input_dumps.append(
                    f"  call dump_array_2d('{arr['name']}.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]})"
                )
            elif arr['dims'] == 3 and len(bounds) == 3:
                input_dumps.append(
                    f"  call dump_array_3d('{arr['name']}.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]}, "
                    f"{bounds[2][0]}, {bounds[2][1]})"
                )
            elif arr['dims'] == 4 and len(bounds) == 4:
                input_dumps.append(
                    f"  call dump_array_4d('{arr['name']}.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]}, "
                    f"{bounds[2][0]}, {bounds[2][1]}, {bounds[3][0]}, {bounds[3][1]})"
                )

    # Dump inout arrays as input
    for arr in args['inout_arrays']:
        if arr['dims'] > 0:
            bounds = extract_array_bounds(arr['dim_str'])
            if arr['dims'] == 2 and len(bounds) == 2:
                input_dumps.append(
                    f"  call dump_array_2d('{arr['name']}_in.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]})"
                )
            elif arr['dims'] == 3 and len(bounds) == 3:
                input_dumps.append(
                    f"  call dump_array_3d('{arr['name']}_in.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]}, "
                    f"{bounds[2][0]}, {bounds[2][1]})"
                )
            elif arr['dims'] == 4 and len(bounds) == 4:
                input_dumps.append(
                    f"  call dump_array_4d('{arr['name']}_in.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]}, "
                    f"{bounds[2][0]}, {bounds[2][1]}, {bounds[3][0]}, {bounds[3][1]})"
                )

    # Dump output and inout arrays as output
    for arr in args['output_arrays'] + args['inout_arrays']:
        if arr['dims'] > 0:
            bounds = extract_array_bounds(arr['dim_str'])
            suffix = '_ref' if arr in args['output_arrays'] else '_ref'
            if arr['dims'] == 2 and len(bounds) == 2:
                output_dumps.append(
                    f"  call dump_array_2d('{arr['name']}_ref.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]})"
                )
            elif arr['dims'] == 3 and len(bounds) == 3:
                output_dumps.append(
                    f"  call dump_array_3d('{arr['name']}_ref.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]}, "
                    f"{bounds[2][0]}, {bounds[2][1]})"
                )
            elif arr['dims'] == 4 and len(bounds) == 4:
                output_dumps.append(
                    f"  call dump_array_4d('{arr['name']}_ref.bin', {arr['name']}, "
                    f"{bounds[0][0]}, {bounds[0][1]}, {bounds[1][0]}, {bounds[1][1]}, "
                    f"{bounds[2][0]}, {bounds[2][1]}, {bounds[3][0]}, {bounds[3][1]})"
                )

    return input_dumps, output_dumps

def add_dump_to_file(filepath, profile_info):
    """Add dump code to a single Fortran file."""

    with open(filepath, 'r') as f:
        content = f.read()

    # Check if already has dump code
    if 'dump_call_count' in content:
        return False, "Already has dump code"

    # Check if has OpenMP parallel
    if '!$omp parallel' not in content.lower():
        return False, "No OpenMP parallel section"

    lines = content.split('\n')
    filename = os.path.basename(filepath)

    # Get profile info for this file
    if filename not in profile_info:
        return False, f"No profile info for {filename}"

    profiles = profile_info[filename]

    # Add use m_dump_kernel to module references
    new_lines = []
    dump_use_added = False

    for i, line in enumerate(lines):
        new_lines.append(line)

        # Add use statement after existing use statements
        if not dump_use_added and re.match(r'\s*use\s+m_comprofile', line, re.IGNORECASE):
            indent = len(line) - len(line.lstrip())
            new_lines.append(' ' * indent + 'use m_dump_kernel')
            dump_use_added = True

    if not dump_use_added:
        # Try to add after any use statement
        for i, line in enumerate(new_lines):
            if re.match(r'\s*use\s+m_\w+', line, re.IGNORECASE):
                indent = len(line) - len(line.lstrip())
                new_lines.insert(i + 1, ' ' * indent + 'use m_dump_kernel')
                dump_use_added = True
                break

    content = '\n'.join(new_lines)
    lines = content.split('\n')

    # Process each OpenMP section
    modified = False
    current_subroutine = None
    section_count = 0

    for profile in profiles:
        subr_name = profile['subroutine']
        target_count = profile['count']
        kernel_name = subr_name.replace('s_', '')

        # Parse arguments for this subroutine
        args = parse_subroutine_args(content, subr_name)
        input_dumps, output_dumps = generate_dump_calls(args, kernel_name)

        # Create unique variable names
        var_suffix = subr_name.replace('s_', '')
        dump_var_decl = f"""
      ! Dump variables
      integer, save :: dump_call_count_{var_suffix} = 0
      integer, parameter :: DUMP_TARGET_{var_suffix} = {target_count}
      logical, save :: dump_done_{var_suffix} = .false.
"""

        # Generate dump code blocks
        input_dump_code = f"""
! Dump input data at target call
dump_call_count_{var_suffix} = dump_call_count_{var_suffix} + 1
if (dump_call_count_{var_suffix} == DUMP_TARGET_{var_suffix} .and. .not. dump_done_{var_suffix}) then
  call dump_init('{kernel_name}')
"""
        for dump_line in input_dumps:
            input_dump_code += dump_line + '\n'
        input_dump_code += "end if\n"

        output_dump_code = f"""
! Dump output data at target call
if (dump_call_count_{var_suffix} == DUMP_TARGET_{var_suffix} .and. .not. dump_done_{var_suffix}) then
"""
        for dump_line in output_dumps:
            output_dump_code += dump_line + '\n'
        output_dump_code += f"""  call dump_finalize()
  dump_done_{var_suffix} = .true.
end if
"""

        # Find and modify the subroutine
        new_lines = []
        in_target_subroutine = False
        var_decl_added = False
        looking_for_omp_parallel = False
        omp_parallel_found = False

        for i, line in enumerate(lines):
            # Track subroutine entry
            subr_match = re.match(rf'\s*subroutine\s+{subr_name}\b', line, re.IGNORECASE)
            if subr_match:
                in_target_subroutine = True
                looking_for_omp_parallel = True

            # Add variable declarations after profiling variables
            if in_target_subroutine and not var_decl_added:
                if 'integer(8) :: loop_len' in line:
                    new_lines.append(line)
                    new_lines.append(dump_var_decl)
                    var_decl_added = True
                    continue

            # Add input dump before !$omp parallel
            if in_target_subroutine and looking_for_omp_parallel and not omp_parallel_found:
                if re.match(r'\s*!\$omp\s+parallel\b', line, re.IGNORECASE):
                    # Add input dump code before this line
                    new_lines.append(input_dump_code)
                    omp_parallel_found = True
                    looking_for_omp_parallel = False

            new_lines.append(line)

            # Add output dump after !$omp end parallel
            if in_target_subroutine and omp_parallel_found:
                if re.match(r'\s*!\$omp\s+end\s+parallel', line, re.IGNORECASE):
                    # Check if next line is profile_stop
                    if i + 1 < len(lines) and 'profile_stop' in lines[i + 1]:
                        # Add after profile_stop
                        new_lines.append(lines[i + 1])
                        new_lines.append(output_dump_code)
                        # Skip the profile_stop line in next iteration
                        lines[i + 1] = '! SKIP_THIS_LINE'
                    else:
                        new_lines.append(output_dump_code)
                    in_target_subroutine = False
                    omp_parallel_found = False
                    modified = True

        lines = [l for l in new_lines if l != '! SKIP_THIS_LINE']
        content = '\n'.join(lines)

    if modified:
        with open(filepath, 'w') as f:
            f.write(content)
        return True, f"Added dump code for {len(profiles)} sections"

    return False, "No modifications made"

def main():
    src_dir = Path('Src')
    profile_path = Path('test_real/omp_profile.txt')

    if not src_dir.exists():
        print("Error: Src directory not found")
        sys.exit(1)

    if not profile_path.exists():
        print("Error: omp_profile.txt not found")
        sys.exit(1)

    # Parse profile data
    print("Parsing omp_profile.txt...")
    profiles = parse_omp_profile(profile_path)
    print(f"Found {len(profiles)} files with profile data")

    # Process each file
    total_modified = 0
    total_skipped = 0

    for filename, file_profiles in sorted(profiles.items()):
        filepath = src_dir / filename
        if not filepath.exists():
            print(f"Warning: {filename} not found")
            continue

        modified, message = add_dump_to_file(filepath, profiles)
        if modified:
            print(f"Modified: {filename} - {message}")
            total_modified += 1
        else:
            print(f"Skipped: {filename} - {message}")
            total_skipped += 1

    print(f"\nTotal: {total_modified} files modified, {total_skipped} files skipped")

if __name__ == '__main__':
    main()
