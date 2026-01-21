#!/usr/bin/env python3
"""
Script to add profiling calls to OpenMP parallel sections in Fortran files.
Adds:
1. Profile ID variables (save attribute)
2. Loop length variable
3. profile_register, profile_start, profile_stop calls
"""

import os
import re
from pathlib import Path

def get_indent(line):
    """Get indentation of a line."""
    return len(line) - len(line.lstrip())

def find_loop_bounds(lines, start_idx, max_search=15):
    """Find loop bounds from do statements following OMP parallel."""
    bounds = []
    for i in range(start_idx, min(start_idx + max_search, len(lines))):
        line = lines[i].strip().lower()
        # Match: do var=start,end
        match = re.match(r'do\s+\w+\s*=\s*([^,]+)\s*,\s*([^\s,!]+)', line)
        if match:
            start, end = match.groups()
            bounds.append((start.strip(), end.strip()))
        # Stop at omp end or other significant construct
        if '!$omp end' in line.lower():
            break
    return bounds

def generate_loop_len_code(bounds, indent):
    """Generate loop length calculation code."""
    if not bounds:
        return f"{indent}loop_len = 1_8\n"

    parts = []
    for start, end in bounds:
        parts.append(f"int(({end})-({start})+1,8)")

    expr = " * ".join(parts)
    if len(expr) > 60:
        # Split into multiple lines
        code = f"{indent}loop_len = "
        for i, part in enumerate(parts):
            if i == 0:
                code += part
            else:
                code += f" &\n{indent}     & * {part}"
        return code + "\n"
    return f"{indent}loop_len = {expr}\n"

def process_file(filepath):
    """Process a single Fortran file."""

    with open(filepath, 'r') as f:
        content = f.read()

    # Check if file has OpenMP parallel sections
    if '!$omp parallel' not in content.lower():
        return False, 0

    # Check if already has profiling
    if 'profile_start' in content:
        return False, 0

    lines = content.split('\n')
    filename = os.path.basename(filepath)

    # Find all subroutines/functions and their OMP sections
    current_subr = None
    subr_start = -1
    subr_end = -1
    sections_info = []

    # First pass: identify subroutines and their OMP sections
    for i, line in enumerate(lines):
        line_lower = line.lower().strip()

        # Track subroutine/function boundaries
        subr_match = re.match(r'\s*subroutine\s+(\w+)', line, re.IGNORECASE)
        func_match = re.match(r'\s*(?:\w+\s+)*function\s+(\w+)', line, re.IGNORECASE)

        if subr_match:
            current_subr = subr_match.group(1)
            subr_start = i
        elif func_match:
            current_subr = func_match.group(1)
            subr_start = i
        elif re.match(r'\s*end\s+subroutine', line, re.IGNORECASE) or \
             re.match(r'\s*end\s+function', line, re.IGNORECASE):
            subr_end = i
            current_subr = None
            subr_start = -1

        # Find OMP parallel sections
        if re.match(r'\s*!\$omp\s+parallel\b', line, re.IGNORECASE):
            if current_subr:
                bounds = find_loop_bounds(lines, i+1)
                sections_info.append({
                    'line_idx': i,
                    'subroutine': current_subr,
                    'bounds': bounds,
                    'indent': get_indent(line)
                })

    if not sections_info:
        return False, 0

    # Group sections by subroutine
    subr_sections = {}
    for sec in sections_info:
        subr = sec['subroutine']
        if subr not in subr_sections:
            subr_sections[subr] = []
        subr_sections[subr].append(sec)

    # Second pass: find where to insert variable declarations
    # Look for the last variable declaration before executable statements
    new_lines = []
    sections_processed = {}  # Track which sections we've processed
    section_counter = {}  # Counter per subroutine

    i = 0
    while i < len(lines):
        line = lines[i]
        line_lower = line.lower().strip()

        # Detect subroutine entry
        subr_match = re.match(r'\s*subroutine\s+(\w+)', line, re.IGNORECASE)
        func_match = re.match(r'\s*(?:\w+\s+)*function\s+(\w+)', line, re.IGNORECASE)

        current_subr_name = None
        if subr_match:
            current_subr_name = subr_match.group(1)
        elif func_match:
            current_subr_name = func_match.group(1)

        if current_subr_name and current_subr_name in subr_sections:
            # We're entering a subroutine that has OMP sections
            # Find the end of variable declarations (look for !-----7 separator)
            new_lines.append(line)
            i += 1

            # Find where to insert profiling variables
            var_insert_idx = -1
            while i < len(lines):
                curr_line = lines[i]
                curr_lower = curr_line.lower().strip()

                new_lines.append(curr_line)

                # Look for the separator line before executable statements
                if re.match(r'!\-+7\-+', curr_lower):
                    var_insert_idx = len(new_lines) - 1
                    i += 1
                    break

                # Also check for common markers of end of declarations
                if ('call ' in curr_lower and not curr_lower.startswith('!')) or \
                   (re.match(r'if\s*\(', curr_lower) and not curr_lower.startswith('!')):
                    var_insert_idx = len(new_lines) - 1
                    break

                i += 1

            # Insert profiling variable declarations
            if var_insert_idx > 0:
                num_sections = len(subr_sections[current_subr_name])
                indent = '      '

                # Build declaration lines
                decl_lines = []
                decl_lines.append(f"\n{indent}! Profiling variables")
                for s_idx in range(1, num_sections + 1):
                    decl_lines.append(f"{indent}integer, save :: prof_id{s_idx} = -1")
                decl_lines.append(f"{indent}integer(8) :: loop_len\n")

                # Insert at var_insert_idx
                new_lines = new_lines[:var_insert_idx] + decl_lines + new_lines[var_insert_idx:]

                # Initialize counter for this subroutine
                section_counter[current_subr_name] = 0

            continue

        # Check for OMP parallel directive
        if re.match(r'\s*!\$omp\s+parallel\b', line, re.IGNORECASE):
            # Find which subroutine this belongs to
            for subr_name, secs in subr_sections.items():
                for sec in secs:
                    if sec['line_idx'] == i and subr_name not in sections_processed:
                        if subr_name not in section_counter:
                            section_counter[subr_name] = 0
                        section_counter[subr_name] += 1
                        sec_num = section_counter[subr_name]

                        indent = ' ' * sec['indent']
                        bounds = sec['bounds']

                        # Add profiling code before OMP parallel
                        prof_code = []
                        prof_code.append(f"\n{indent}! Register profiling section (first call only)")
                        prof_code.append(f"{indent}if (prof_id{sec_num} < 0) then")
                        prof_code.append(f"{indent}  prof_id{sec_num} = profile_register('{filename}', '{subr_name}', &")
                        prof_code.append(f"{indent}   & 'OMP section {sec_num}')")
                        prof_code.append(f"{indent}end if")
                        prof_code.append(generate_loop_len_code(bounds, indent).rstrip())
                        prof_code.append(f"{indent}call profile_start(prof_id{sec_num})")
                        prof_code.append("")

                        for pc_line in prof_code:
                            new_lines.append(pc_line)

                        new_lines.append(line)
                        i += 1

                        # Now find the matching !$omp end parallel and add profile_stop
                        omp_depth = 1
                        while i < len(lines) and omp_depth > 0:
                            curr_line = lines[i]
                            curr_lower = curr_line.lower().strip()

                            if re.match(r'!\$omp\s+parallel\b', curr_lower):
                                omp_depth += 1
                            elif re.match(r'!\$omp\s+end\s+parallel\b', curr_lower):
                                omp_depth -= 1

                            new_lines.append(curr_line)
                            i += 1

                            if omp_depth == 0:
                                # Add profile_stop after end parallel
                                new_lines.append(f"\n{indent}call profile_stop(prof_id{sec_num}, loop_len)")
                                break

                        if subr_name in sections_processed:
                            sections_processed[subr_name] += 1
                        else:
                            sections_processed[subr_name] = 1
                        break
                else:
                    continue
                break
            else:
                new_lines.append(line)
                i += 1
        else:
            new_lines.append(line)
            i += 1

    with open(filepath, 'w') as f:
        f.write('\n'.join(new_lines))

    return True, len(sections_info)

def main():
    src_dir = Path('Src')

    if not src_dir.exists():
        print("Error: Src directory not found")
        return

    total_files = 0
    total_sections = 0

    for f90_file in sorted(src_dir.glob('*.f90')):
        try:
            modified, sections = process_file(f90_file)
            if modified:
                print(f"Processed: {f90_file.name} ({sections} sections)")
                total_files += 1
                total_sections += sections
        except Exception as e:
            print(f"Error processing {f90_file.name}: {e}")

    print(f"\nTotal: {total_files} files, {total_sections} sections")

if __name__ == '__main__':
    main()
