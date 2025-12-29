#!/usr/bin/env python3
"""
Script to add 'use m_comprofile' to all Fortran files with OpenMP parallel sections.
"""

import os
import re
from pathlib import Path

def add_use_statement(filepath):
    """Add use m_comprofile to a Fortran file."""

    with open(filepath, 'r') as f:
        content = f.read()

    # Check if file has OpenMP parallel sections
    if '!$omp parallel' not in content.lower():
        return False

    # Check if already has profiling module
    if 'use m_comprofile' in content.lower():
        return False

    lines = content.split('\n')
    new_lines = []
    use_added = False

    for i, line in enumerate(lines):
        new_lines.append(line)

        # Find first 'use m_' statement and add after it
        if not use_added and re.match(r'\s*use\s+m_\w+', line, re.IGNORECASE):
            # Get indentation from current line
            indent = len(line) - len(line.lstrip())
            new_lines.append(' ' * indent + 'use m_comprofile')
            use_added = True

    if not use_added:
        return False

    with open(filepath, 'w') as f:
        f.write('\n'.join(new_lines))

    return True

def main():
    src_dir = Path('Src')

    if not src_dir.exists():
        print("Error: Src directory not found")
        return

    modified_count = 0

    for f90_file in sorted(src_dir.glob('*.f90')):
        if add_use_statement(f90_file):
            print(f"Added use m_comprofile to: {f90_file.name}")
            modified_count += 1

    print(f"\nTotal files modified: {modified_count}")

if __name__ == '__main__':
    main()
