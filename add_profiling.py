#!/usr/bin/env python3
"""
Script to add profiling code to OpenMP parallel sections in Fortran files.
"""

import os
import re
import sys
from pathlib import Path

def parse_loop_bounds(content, omp_start_line):
    """Extract loop bounds from do loops following OpenMP parallel directive."""
    lines = content.split('\n')
    loop_info = []

    # Look for do loops in the next 20 lines
    for i in range(omp_start_line, min(omp_start_line + 20, len(lines))):
        line = lines[i].strip().lower()
        # Match: do var=start,end
        match = re.match(r'do\s+(\w+)\s*=\s*([^,]+)\s*,\s*([^\s]+)', line)
        if match:
            var, start, end = match.groups()
            loop_info.append((var.strip(), start.strip(), end.strip()))
        elif line.startswith('!$omp end'):
            break

    return loop_info

def generate_loop_len_expr(loop_info):
    """Generate expression to calculate loop length."""
    if not loop_info:
        return "1_8"

    parts = []
    for var, start, end in loop_info:
        # Generate expression: (end - start + 1)
        parts.append(f"int(({end})-({start})+1,8)")

    return " * ".join(parts)

def add_profiling_to_file(filepath):
    """Add profiling code to a single Fortran file."""

    with open(filepath, 'r') as f:
        content = f.read()

    # Check if file has OpenMP parallel sections
    if '!$omp parallel' not in content.lower():
        return False, 0

    # Check if already has profiling
    if 'use m_comprofile' in content:
        return False, 0

    lines = content.split('\n')
    filename = os.path.basename(filepath)

    # Find module/subroutine/function names
    current_subr = "unknown"

    # Find all OpenMP parallel sections
    omp_sections = []
    for i, line in enumerate(lines):
        # Track current subroutine
        subr_match = re.match(r'\s*subroutine\s+(\w+)', line, re.IGNORECASE)
        if subr_match:
            current_subr = subr_match.group(1)
        func_match = re.match(r'\s*function\s+(\w+)', line, re.IGNORECASE)
        if func_match:
            current_subr = func_match.group(1)

        # Find OpenMP parallel directive
        if re.match(r'\s*!\$omp\s+parallel\b', line, re.IGNORECASE):
            omp_sections.append({
                'line_num': i,
                'subroutine': current_subr,
                'line': line
            })

    if not omp_sections:
        return False, 0

    # Add use m_comprofile to module references
    new_lines = []
    use_added = False

    for i, line in enumerate(lines):
        new_lines.append(line)

        # Add use statement after first "use" statement
        if not use_added and re.match(r'\s*use\s+m_\w+', line, re.IGNORECASE):
            # Check if next line is also a use statement
            if i + 1 < len(lines) and not re.match(r'\s*use\s+', lines[i+1], re.IGNORECASE):
                # This is the last use statement, add after it
                pass
            else:
                continue
            # Get indentation
            indent = len(line) - len(line.lstrip())
            new_lines.insert(-1, ' ' * indent + 'use m_comprofile')
            use_added = True

    # If no use statement was found, try adding after 'use m_' pattern
    if not use_added:
        for i, line in enumerate(new_lines):
            if re.match(r'\s*use\s+m_', line, re.IGNORECASE):
                indent = len(line) - len(line.lstrip())
                new_lines.insert(i, ' ' * indent + 'use m_comprofile')
                use_added = True
                break

    content = '\n'.join(new_lines)

    # Now process each OpenMP section (work backwards to preserve line numbers)
    section_count = len(omp_sections)

    # For each subroutine, we need to add variable declarations
    subroutines_processed = set()

    for idx, section in enumerate(reversed(omp_sections)):
        section_id = section_count - idx
        subr = section['subroutine']

        # Create unique profile ID variable name
        prof_var = f"prof_id{section_id}"

    # Write the modified content
    with open(filepath, 'w') as f:
        f.write(content)

    return True, len(omp_sections)

def main():
    src_dir = Path('Src')

    if not src_dir.exists():
        print("Error: Src directory not found")
        sys.exit(1)

    total_files = 0
    total_sections = 0

    for f90_file in sorted(src_dir.glob('*.f90')):
        modified, sections = add_profiling_to_file(f90_file)
        if modified:
            print(f"Modified: {f90_file.name} ({sections} sections)")
            total_files += 1
            total_sections += sections

    print(f"\nTotal: {total_files} files, {total_sections} sections")

if __name__ == '__main__':
    main()
