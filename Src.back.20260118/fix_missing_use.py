#!/usr/bin/env python3
"""
Fix missing 'use m_comprofile' in files that have profile calls but no use statement.
"""

import re
from pathlib import Path

def fix_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    # Check if has profile calls but no use statement
    if 'use m_comprofile' in content.lower():
        return False
    if 'profile_start' not in content and 'profile_register' not in content:
        return False

    lines = content.split('\n')
    new_lines = []
    added = False

    for i, line in enumerate(lines):
        new_lines.append(line)

        # Look for "! Module reference" section
        if not added and '! Module reference' in line:
            # Check next few lines for 'none' or existing use
            for j in range(i+1, min(i+5, len(lines))):
                if lines[j].strip().lower() == '!     none':
                    # Replace 'none' with use statement
                    # Find position and insert
                    pass
                elif 'use ' in lines[j].lower():
                    break

        # Add after implicit none if we haven't added yet
        if not added and re.match(r'\s*implicit\s+none', line, re.IGNORECASE):
            # Find indentation
            indent = '      '
            # Look back for module reference section
            for j in range(len(new_lines)-1, max(0, len(new_lines)-20), -1):
                if '! Module reference' in new_lines[j]:
                    # Insert use statement after Module reference comment
                    for k in range(j+1, len(new_lines)):
                        if new_lines[k].strip().lower() == '!     none':
                            new_lines[k] = f'{indent}use m_comprofile'
                            added = True
                            break
                        elif new_lines[k].strip() == '' or new_lines[k].strip().startswith('!'):
                            continue
                        else:
                            break
                    break

    if not added:
        # Alternative: add before implicit none
        new_lines2 = []
        for i, line in enumerate(lines):
            if not added and re.match(r'\s*implicit\s+none', line, re.IGNORECASE):
                indent = '      '
                new_lines2.append(f'{indent}use m_comprofile')
                new_lines2.append('')
                added = True
            new_lines2.append(line)
        new_lines = new_lines2

    if added:
        with open(filepath, 'w') as f:
            f.write('\n'.join(new_lines))
        return True
    return False

def main():
    count = 0
    for f90_file in sorted(Path('.').glob('*.f90')):
        if fix_file(f90_file):
            print(f"Fixed: {f90_file.name}")
            count += 1
    print(f"\nTotal fixed: {count}")

if __name__ == '__main__':
    main()
