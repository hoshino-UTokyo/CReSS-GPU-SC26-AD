#!/usr/bin/env python3
"""Update README.md files with Dump Data Requirements section."""

import os
import re
import glob

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(BASE_DIR, 'Src')
KERNEL_DIR = os.path.join(BASE_DIR, 'Kernel_benchmark')


def extract_dump_info(srcfile):
    """Extract dump information from a source file."""
    if not os.path.exists(srcfile):
        return None

    with open(srcfile, 'r') as f:
        content = f.read()

    if 'call dump_init' not in content:
        return None

    info = {
        'kernel_name': None,
        'target_call': None,
        'scalars': [],
        'input_arrays': [],
        'output_arrays': []
    }

    # Extract kernel name
    match = re.search(r"call dump_init\('([^']+)'\)", content)
    if match:
        info['kernel_name'] = match.group(1)

    # Extract target call count
    match = re.search(r'DUMP_TARGET_\w+\s*=\s*(\d+)', content)
    if match:
        info['target_call'] = match.group(1)

    # Extract scalar dumps
    for match in re.finditer(r"call dump_scalar_([irc])\('([^']+)'", content):
        scalar_type = {'i': 'integer', 'r': 'real', 'c': 'character'}[match.group(1)]
        info['scalars'].append((match.group(2), scalar_type))

    # Extract array dumps
    for match in re.finditer(r"call dump_array_(\d+)d\('([^']+\.bin)',\s*([^,]+),\s*([^)]+)\)", content):
        dims = match.group(1)
        filename = match.group(2)
        varname = match.group(3).strip()
        bounds = match.group(4).strip()

        if '_ref.' in filename:
            info['output_arrays'].append((filename, varname, dims, bounds))
        else:
            info['input_arrays'].append((filename, varname, dims, bounds))

    return info


def generate_dump_section(info):
    """Generate markdown section for dump data requirements."""
    if not info or not info['kernel_name']:
        return None

    lines = []
    lines.append('')
    lines.append('## Dump Data Requirements')
    lines.append('')
    lines.append(f"**Kernel Name**: `{info['kernel_name']}`")
    lines.append(f"**Dump Target Call**: {info['target_call']}")
    lines.append(f"**Dump Directory**: `test_real/kernel_dump/{info['kernel_name']}/`")
    lines.append('')

    if info['scalars']:
        lines.append('### Input Scalars (params.txt)')
        lines.append('| Name | Type |')
        lines.append('|------|------|')
        for name, stype in info['scalars']:
            lines.append(f'| `{name}` | {stype} |')
        lines.append('')

    if info['input_arrays']:
        lines.append('### Input Arrays')
        lines.append('| File | Variable | Dims | Bounds |')
        lines.append('|------|----------|------|--------|')
        for filename, varname, dims, bounds in info['input_arrays']:
            lines.append(f'| `{filename}` | `{varname}` | {dims}D | `{bounds}` |')
        lines.append('')

    if info['output_arrays']:
        lines.append('### Output Arrays (Reference)')
        lines.append('| File | Variable | Dims | Bounds |')
        lines.append('|------|----------|------|--------|')
        for filename, varname, dims, bounds in info['output_arrays']:
            lines.append(f'| `{filename}` | `{varname}` | {dims}D | `{bounds}` |')
        lines.append('')

    return '\n'.join(lines)


def update_readme(readme_path, dump_section):
    """Update README.md with dump section."""
    with open(readme_path, 'r') as f:
        content = f.read()

    # Check if already has dump section
    if '## Dump Data Requirements' in content:
        return False

    # Append dump section
    with open(readme_path, 'a') as f:
        f.write(dump_section)

    return True


def main():
    updated = 0
    skipped = 0

    for readme in glob.glob(os.path.join(KERNEL_DIR, '*/README.md')):
        kernel_dir = os.path.dirname(readme)

        # Read README to find source file
        with open(readme, 'r') as f:
            readme_content = f.read()

        match = re.search(r'Src/([a-z0-9_]+\.f90)', readme_content)
        if not match:
            continue

        srcfile = os.path.join(SRC_DIR, match.group(1))
        info = extract_dump_info(srcfile)

        if not info:
            continue

        dump_section = generate_dump_section(info)
        if not dump_section:
            continue

        if update_readme(readme, dump_section):
            print(f"Updated: {readme}")
            updated += 1
        else:
            skipped += 1

    print(f"\nTotal updated: {updated}, Already had dump section: {skipped}")


if __name__ == '__main__':
    main()
