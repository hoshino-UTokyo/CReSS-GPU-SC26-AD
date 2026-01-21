#!/usr/bin/env python3
"""
Add runtime profiling information to meta_info annotations.
Reads omp_profile.txt and updates source files with runtime stats.
"""

import re
import os
from collections import defaultdict

def parse_profile(profile_path):
    """Parse omp_profile.txt and return dict of profiling data."""
    profiles = {}
    total_time = 0.0

    with open(profile_path, 'r') as f:
        lines = f.readlines()

    # Skip header lines
    in_data = False
    for line in lines:
        line = line.strip()

        # Skip header and separator lines
        if line.startswith('===') or line.startswith('---'):
            continue
        if 'ID' in line and 'File' in line and 'Subroutine' in line:
            in_data = True
            continue
        if 'Unexecuted' in line:
            break  # Stop at unexecuted sections

        if in_data and line:
            # Parse: ID  File  Subroutine  Count  AvgLoops  TotalTime  AvgTime
            parts = line.split()
            if len(parts) >= 7:
                try:
                    section_id = int(parts[0])
                    filename = parts[1]
                    subroutine = parts[2]
                    count = int(parts[3])
                    avg_loops = float(parts[4])
                    total_time_sec = float(parts[5])
                    avg_time_ms = float(parts[6])

                    # Key by filename (without path)
                    key = filename.lower()
                    if key not in profiles:
                        profiles[key] = []

                    profiles[key].append({
                        'id': section_id,
                        'filename': filename,
                        'subroutine': subroutine,
                        'count': count,
                        'avg_loops': avg_loops,
                        'total_time': total_time_sec,
                        'avg_time': avg_time_ms
                    })

                    total_time += total_time_sec
                except (ValueError, IndexError):
                    pass

    return profiles, total_time


def format_runtime_info(profile_data, total_time):
    """Format runtime info as annotation text."""
    lines = []
    lines.append("! Runtime:")

    pct = (profile_data['total_time'] / total_time * 100) if total_time > 0 else 0

    lines.append(f"!   - Calls: {profile_data['count']}")

    # Format avg_loops nicely
    avg_loops = profile_data['avg_loops']
    if avg_loops >= 1e6:
        loops_str = f"{avg_loops/1e6:.1f}M"
    elif avg_loops >= 1e3:
        loops_str = f"{avg_loops/1e3:.1f}K"
    else:
        loops_str = f"{avg_loops:.0f}"
    lines.append(f"!   - AvgLoops: {loops_str}")

    lines.append(f"!   - TotalTime: {profile_data['total_time']:.3f}s ({pct:.2f}%)")
    lines.append(f"!   - AvgTime: {profile_data['avg_time']:.3f}ms")

    return '\n'.join(lines)


def update_source_file(filepath, profiles, total_time):
    """Update meta_info annotations in a source file with runtime info."""

    filename = os.path.basename(filepath).lower()

    if filename not in profiles:
        return False, "No profile data for this file"

    profile_list = profiles[filename]

    with open(filepath, 'r') as f:
        content = f.read()

    # Check if file has meta_info annotations
    if '!@llm start meta_info' not in content:
        return False, "No meta_info annotation found"

    # Check if Runtime info already exists
    if '! Runtime:' in content:
        return False, "Runtime info already exists"

    lines = content.split('\n')
    new_lines = []
    profile_idx = 0
    modified = False

    i = 0
    while i < len(lines):
        line = lines[i]
        new_lines.append(line)

        # Look for end of meta_info block to insert runtime info before it
        if '!@llm end meta_info' in line:
            # Insert runtime info before the end marker
            if profile_idx < len(profile_list):
                runtime_info = format_runtime_info(profile_list[profile_idx], total_time)
                # Insert before the end marker
                new_lines.pop()  # Remove the end marker line
                new_lines.append(runtime_info)
                new_lines.append(line)  # Re-add the end marker
                profile_idx += 1
                modified = True

        i += 1

    if modified:
        with open(filepath, 'w') as f:
            f.write('\n'.join(new_lines))
        return True, f"Updated {profile_idx} annotation(s)"

    return False, "No modifications made"


def main():
    profile_path = 'test_real/omp_profile.txt'
    src_dir = 'Src'

    print("Parsing profile data...")
    profiles, total_time = parse_profile(profile_path)
    print(f"Found {sum(len(v) for v in profiles.values())} profiled sections")
    print(f"Total execution time: {total_time:.3f}s")
    print()

    # Find all Fortran files
    updated_count = 0
    skipped_count = 0

    for filename in sorted(os.listdir(src_dir)):
        if filename.endswith('.f90') or filename.endswith('.F90'):
            filepath = os.path.join(src_dir, filename)
            success, msg = update_source_file(filepath, profiles, total_time)
            if success:
                print(f"[OK] {filename}: {msg}")
                updated_count += 1
            else:
                skipped_count += 1

    print()
    print(f"Updated: {updated_count} files")
    print(f"Skipped: {skipped_count} files")


if __name__ == '__main__':
    main()
