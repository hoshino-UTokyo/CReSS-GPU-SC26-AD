#!/bin/bash
#
# setup_data_symlinks.sh
#
# Replace data directories with symbolic links to test_real/kernel_dump/
# This avoids copying dump data and ensures benchmarks use the latest data.
#
# Usage:
#   ./setup_data_symlinks.sh          # Create/update symlinks
#   ./setup_data_symlinks.sh --check  # Check symlink status only
#

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
KERNEL_DUMP_DIR="../../test_real/kernel_dump"
KERNEL_DUMP_ABS="$SCRIPT_DIR/../test_real/kernel_dump"

check_only=false
if [ "$1" = "--check" ]; then
    check_only=true
fi

cd "$SCRIPT_DIR"

created=0
skipped=0
missing=0

for bench_dir in */; do
    bench_dir=${bench_dir%/}

    # Skip non-benchmark directories
    [[ ! "$bench_dir" =~ ^[0-9] ]] && continue

    # Extract kernel name from directory name
    # Patterns: 012_advbspi_subroutine, 087_disptke_s_disptke, 009_adjstuv_subroutine_sec1
    kernel_name=$(echo "$bench_dir" | sed 's/^[0-9]*_//' | sed 's/_subroutine.*$//' | sed 's/_s_.*$//')

    # Check if kernel dump exists
    if [ ! -d "$KERNEL_DUMP_ABS/$kernel_name" ]; then
        ((missing++))
        continue
    fi

    data_path="$bench_dir/data"

    if $check_only; then
        if [ -L "$data_path" ]; then
            target=$(readlink "$data_path")
            echo "[OK] $bench_dir -> $target"
            ((created++))
        elif [ -d "$data_path" ]; then
            echo "[DIR] $bench_dir (not a symlink)"
            ((skipped++))
        else
            echo "[MISSING] $bench_dir"
            ((missing++))
        fi
    else
        # Remove existing data (file, directory, or symlink)
        if [ -e "$data_path" ] || [ -L "$data_path" ]; then
            rm -rf "$data_path"
        fi

        # Create symlink
        ln -s "$KERNEL_DUMP_DIR/$kernel_name" "$data_path"
        echo "[$bench_dir] -> kernel_dump/$kernel_name"
        ((created++))
    fi
done

echo ""
echo "=== Summary ==="
if $check_only; then
    echo "Symlinks:    $created"
    echo "Directories: $skipped"
    echo "Missing:     $missing"
else
    echo "Created: $created symlinks"
    echo "Missing kernel dumps: $missing (no symlink created)"
fi
