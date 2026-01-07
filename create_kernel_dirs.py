#!/usr/bin/env python3
"""
Create Kernel_benchmark directory structure for all 387 OpenMP parallel sections.
Extracts annotations directly from source files.
"""

import os
import re
import glob

def parse_profile(profile_path):
    """Parse omp_profile.txt and return dict keyed by (filename, subroutine)."""
    profiles = {}
    total_time = 0.0

    with open(profile_path, 'r') as f:
        lines = f.readlines()

    in_data = False
    in_unexecuted = False

    for line in lines:
        line_stripped = line.strip()

        if line_stripped.startswith('===') or line_stripped.startswith('---'):
            continue
        if 'ID' in line_stripped and 'File' in line_stripped:
            in_data = True
            continue
        if 'Unexecuted' in line_stripped:
            in_unexecuted = True
            continue

        if in_data and line_stripped and '.f90' in line_stripped.lower():
            parts = line_stripped.split()
            if len(parts) >= 7:
                try:
                    section_id = int(parts[0])
                    filename = parts[1].lower()
                    subroutine = parts[2].lower()
                    count = int(parts[3])
                    avg_loops = float(parts[4])
                    total_time_sec = float(parts[5])
                    avg_time = float(parts[6])

                    key = (filename, subroutine)
                    if key not in profiles:
                        profiles[key] = []

                    profiles[key].append({
                        'id': section_id,
                        'count': count,
                        'avg_loops': avg_loops,
                        'total_time': total_time_sec,
                        'avg_time': avg_time,
                        'executed': not in_unexecuted
                    })

                    if not in_unexecuted:
                        total_time += total_time_sec
                except (ValueError, IndexError):
                    pass

    return profiles, total_time


def extract_annotations(src_dir):
    """Extract all meta_info annotations from source files."""
    annotations = []

    for filepath in sorted(glob.glob(os.path.join(src_dir, '*.f90'))):
        filename = os.path.basename(filepath)

        with open(filepath, 'r') as f:
            content = f.read()

        # Find all meta_info blocks
        pattern = r'!@llm start meta_info.*?!@llm end meta_info[^\n]*'
        matches = list(re.finditer(pattern, content, re.DOTALL))

        for idx, match in enumerate(matches):
            block = match.group(0)
            start_pos = match.start()

            # Extract location info
            loc_match = re.search(r'Location:\s*(\S+)\s*::\s*(\S+)', block)
            if loc_match:
                subroutine = loc_match.group(2)
            else:
                subroutine = 'unknown'

            # Extract summary
            sum_match = re.search(r'Summary\s*:\s*(.+?)(?=\n!)', block, re.DOTALL)
            summary = sum_match.group(1).strip() if sum_match else ''
            # Clean up multiline summary
            summary = ' '.join(summary.replace('\n', ' ').replace('!', '').split())

            # Extract GPU difficulty
            diff_match = re.search(r'GPU diff:\s*(\S+)', block)
            difficulty = diff_match.group(1) if diff_match else 'Unknown'

            # Extract Runtime info if present
            runtime_match = re.search(r'Runtime:.*?Calls:\s*(\d+).*?AvgLoops:\s*(\S+).*?TotalTime:\s*(\S+).*?AvgTime:\s*(\S+)', block, re.DOTALL)
            runtime = None
            if runtime_match:
                runtime = {
                    'calls': int(runtime_match.group(1)),
                    'avg_loops': runtime_match.group(2),
                    'total_time': runtime_match.group(3),
                    'avg_time': runtime_match.group(4)
                }

            # Calculate line number
            line_number = content[:start_pos].count('\n') + 1

            annotations.append({
                'filename': filename,
                'subroutine': subroutine,
                'section_idx': idx + 1,  # 1-based section index within file
                'total_sections': len(matches),
                'line_number': line_number,
                'summary': summary[:200],  # Truncate long summaries
                'difficulty': difficulty,
                'runtime': runtime,
                'block': block
            })

    return annotations


def create_readme(ann, kernel_id, total_time):
    """Create README content for a kernel directory."""
    lines = []
    lines.append(f"# Kernel {kernel_id:03d}: {ann['subroutine']}")
    lines.append("")
    lines.append("## Source Location")
    lines.append(f"- **File**: Src/{ann['filename']}")
    lines.append(f"- **Subroutine**: {ann['subroutine']}")
    lines.append(f"- **Line**: ~{ann['line_number']}")
    if ann['total_sections'] > 1:
        lines.append(f"- **Section**: {ann['section_idx']} of {ann['total_sections']} in this subroutine")
    lines.append("")

    lines.append("## Analysis")
    lines.append(f"- **GPU Difficulty**: {ann['difficulty']}")
    if ann['summary']:
        lines.append(f"- **Summary**: {ann['summary']}")
    lines.append("")

    if ann['runtime']:
        rt = ann['runtime']
        lines.append("## Runtime Profile (from test_real)")
        lines.append(f"- **Calls**: {rt['calls']}")
        lines.append(f"- **Average Loop Length**: {rt['avg_loops']}")
        lines.append(f"- **Total Time**: {rt['total_time']}")
        lines.append(f"- **Average Time per Call**: {rt['avg_time']}")
    else:
        lines.append("## Runtime Profile")
        lines.append("- **Status**: Not executed or no runtime data available")
    lines.append("")

    lines.append("## Original Annotation")
    lines.append("```fortran")
    for line in ann['block'].split('\n'):
        lines.append(line)
    lines.append("```")
    lines.append("")

    lines.append("## Benchmark Files (TODO)")
    lines.append("- `kernel.f90` - Extracted kernel code")
    lines.append("- `driver.f90` - Benchmark driver")
    lines.append("- `data/` - Input data for benchmark")
    lines.append("- `Makefile` - Build configuration")
    lines.append("")

    return '\n'.join(lines)


def main():
    profile_path = 'test_real/omp_profile.txt'
    src_dir = 'Src'
    base_dir = 'Kernel_benchmark'

    print("Parsing profile data...")
    profiles, total_time = parse_profile(profile_path)
    print(f"Total execution time: {total_time:.3f}s")

    print("Extracting annotations from source files...")
    annotations = extract_annotations(src_dir)
    print(f"Found {len(annotations)} annotations")

    # Create base directory
    os.makedirs(base_dir, exist_ok=True)

    # Statistics
    difficulty_counts = {'Easy': 0, 'Medium': 0, 'Hard': 0, 'Unknown': 0}
    executed_count = 0

    # Create directories
    for kernel_id, ann in enumerate(annotations, 1):
        # Directory name format: {ID}_{filename}_{subroutine}[_sec{N}]
        filename_base = ann['filename'].replace('.f90', '').replace('.F90', '')
        if ann['total_sections'] > 1:
            dir_name = f"{kernel_id:03d}_{filename_base}_{ann['subroutine']}_sec{ann['section_idx']}"
        else:
            dir_name = f"{kernel_id:03d}_{filename_base}_{ann['subroutine']}"

        dir_path = os.path.join(base_dir, dir_name)
        os.makedirs(dir_path, exist_ok=True)

        # Create README
        readme_content = create_readme(ann, kernel_id, total_time)
        readme_path = os.path.join(dir_path, 'README.md')
        with open(readme_path, 'w') as f:
            f.write(readme_content)

        # Update statistics
        diff = ann['difficulty']
        if diff in difficulty_counts:
            difficulty_counts[diff] += 1
        else:
            difficulty_counts['Unknown'] += 1

        if ann['runtime']:
            executed_count += 1

    # Create main index
    index_lines = []
    index_lines.append("# Kernel Benchmark Directory")
    index_lines.append("")
    index_lines.append("Individual benchmark setups for each OpenMP parallel section.")
    index_lines.append("")
    index_lines.append("## Statistics")
    index_lines.append(f"- **Total Kernels**: {len(annotations)}")
    index_lines.append(f"- **With Runtime Data**: {executed_count}")
    index_lines.append(f"- **Total Execution Time**: {total_time:.3f}s")
    index_lines.append("")
    index_lines.append("## GPU Difficulty Distribution")
    index_lines.append(f"- Easy: {difficulty_counts['Easy']}")
    index_lines.append(f"- Medium: {difficulty_counts['Medium']}")
    index_lines.append(f"- Hard: {difficulty_counts['Hard']}")
    index_lines.append("")
    index_lines.append("## Directory Naming Convention")
    index_lines.append("```")
    index_lines.append("{ID:03d}_{filename}_{subroutine}[_sec{N}]/")
    index_lines.append("```")
    index_lines.append("- `_sec{N}` suffix added when subroutine has multiple OpenMP sections")
    index_lines.append("")

    # Top kernels by difficulty
    index_lines.append("## Kernels by Difficulty")
    index_lines.append("")
    for diff in ['Easy', 'Medium', 'Hard']:
        index_lines.append(f"### {diff} ({difficulty_counts[diff]})")
        index_lines.append("")
        count = 0
        for kernel_id, ann in enumerate(annotations, 1):
            if ann['difficulty'] == diff:
                filename_base = ann['filename'].replace('.f90', '')
                if ann['total_sections'] > 1:
                    dir_name = f"{kernel_id:03d}_{filename_base}_{ann['subroutine']}_sec{ann['section_idx']}"
                else:
                    dir_name = f"{kernel_id:03d}_{filename_base}_{ann['subroutine']}"
                index_lines.append(f"- [{dir_name}](./{dir_name}/)")
                count += 1
                if count >= 50 and diff == 'Easy':
                    index_lines.append(f"- ... and {difficulty_counts[diff] - 50} more")
                    break
        index_lines.append("")

    index_path = os.path.join(base_dir, 'README.md')
    with open(index_path, 'w') as f:
        f.write('\n'.join(index_lines))

    print(f"\nCreated {len(annotations)} directories in {base_dir}/")
    print(f"  Easy: {difficulty_counts['Easy']}")
    print(f"  Medium: {difficulty_counts['Medium']}")
    print(f"  Hard: {difficulty_counts['Hard']}")
    print(f"\nIndex created: {index_path}")


if __name__ == '__main__':
    main()
