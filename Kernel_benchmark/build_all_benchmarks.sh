#!/bin/bash
#==============================================================================
# Build All Kernel Benchmarks
#==============================================================================
# Usage: ./build_all_benchmarks.sh [OPTIONS]
#
# Options:
#   -c          Clean before building (make clean && make)
#   -j N        Parallel build jobs (default: 1)
#   -q          Quiet mode (only show summary)
#   -s          Stop on first build failure
#   -f FILTER   Only build kernels matching pattern (e.g., "turb*")
#   -h          Show this help
#==============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Default settings
CLEAN=0
JOBS=1
QUIET=0
STOP_ON_FAIL=0
FILTER="*"

# Parse arguments
while getopts "cj:qsf:h" opt; do
    case $opt in
        c) CLEAN=1 ;;
        j) JOBS=$OPTARG ;;
        q) QUIET=1 ;;
        s) STOP_ON_FAIL=1 ;;
        f) FILTER=$OPTARG ;;
        h)
            head -17 "$0" | tail -16
            exit 0
            ;;
        *)
            echo "Unknown option: -$opt"
            exit 1
            ;;
    esac
done

# Statistics
TOTAL=0
BUILT=0
FAILED=0
SKIPPED=0
NO_MAKEFILE=0
NO_SOURCE=0
NO_DATA=0

# Results storage
BUILT_LIST=""
FAILED_LIST=""
SKIPPED_LIST=""

# Header
echo "=============================================================================="
echo " Kernel Benchmark Builder"
echo "=============================================================================="
echo " Directory: $SCRIPT_DIR"
echo " Clean build: $([ $CLEAN -eq 1 ] && echo 'Yes' || echo 'No')"
echo " Filter: $FILTER"
echo " Date: $(date)"
echo "=============================================================================="
echo ""

# Find all kernel directories matching the filter
for dir in [0-9][0-9][0-9]_*/; do
    [ -d "$dir" ] || continue

    kernel_name=$(basename "$dir")

    # Apply filter
    if [[ ! "$kernel_name" == $FILTER ]]; then
        continue
    fi

    TOTAL=$((TOTAL + 1))

    # Check if Makefile exists
    if [ ! -f "$dir/Makefile" ]; then
        if [ $QUIET -eq 0 ]; then
            echo "[$kernel_name] SKIP: No Makefile"
        fi
        NO_MAKEFILE=$((NO_MAKEFILE + 1))
        SKIPPED=$((SKIPPED + 1))
        SKIPPED_LIST="$SKIPPED_LIST $kernel_name(no-makefile)"
        continue
    fi

    # Check if source file exists
    if [ ! -f "$dir/kernel_benchmark.f90" ]; then
        if [ $QUIET -eq 0 ]; then
            echo "[$kernel_name] SKIP: No kernel_benchmark.f90"
        fi
        NO_SOURCE=$((NO_SOURCE + 1))
        SKIPPED=$((SKIPPED + 1))
        SKIPPED_LIST="$SKIPPED_LIST $kernel_name(no-source)"
        continue
    fi

    # Check if data directory exists (symlink or directory)
    if [ ! -d "$dir/data" ] && [ ! -L "$dir/data" ]; then
        if [ $QUIET -eq 0 ]; then
            echo "[$kernel_name] SKIP: No data directory (symlink or directory)"
        fi
        NO_DATA=$((NO_DATA + 1))
        SKIPPED=$((SKIPPED + 1))
        SKIPPED_LIST="$SKIPPED_LIST $kernel_name(no-data)"
        continue
    fi

    # Check if data symlink target exists
    if [ -L "$dir/data" ] && [ ! -e "$dir/data" ]; then
        if [ $QUIET -eq 0 ]; then
            echo "[$kernel_name] SKIP: Data symlink broken"
        fi
        NO_DATA=$((NO_DATA + 1))
        SKIPPED=$((SKIPPED + 1))
        SKIPPED_LIST="$SKIPPED_LIST $kernel_name(broken-symlink)"
        continue
    fi

    # Build
    if [ $QUIET -eq 0 ]; then
        echo -n "[$kernel_name] Building... "
    fi

    cd "$dir"

    # Run make
    if [ $CLEAN -eq 1 ]; then
        make_output=$(make clean 2>&1 && make 2>&1)
    else
        make_output=$(make 2>&1)
    fi
    make_exit=$?

    cd "$SCRIPT_DIR"

    if [ $make_exit -eq 0 ]; then
        BUILT=$((BUILT + 1))
        BUILT_LIST="$BUILT_LIST $kernel_name"
        if [ $QUIET -eq 0 ]; then
            echo "OK"
        fi
    else
        FAILED=$((FAILED + 1))
        FAILED_LIST="$FAILED_LIST $kernel_name"
        if [ $QUIET -eq 0 ]; then
            echo "FAILED"
            echo "--- Build output ---"
            echo "$make_output"
            echo "--------------------"
        fi

        if [ $STOP_ON_FAIL -eq 1 ]; then
            echo ""
            echo "Stopping on first build failure."
            break
        fi
    fi
done

# Summary
echo ""
echo "=============================================================================="
echo " Build Summary"
echo "=============================================================================="
echo " Total directories:     $TOTAL"
echo " Successfully built:    $BUILT"
echo " Build failed:          $FAILED"
echo " Skipped:               $SKIPPED"
echo "   - No Makefile:       $NO_MAKEFILE"
echo "   - No source file:    $NO_SOURCE"
echo "   - No data:           $NO_DATA"
echo "=============================================================================="

if [ -n "$BUILT_LIST" ] && [ $QUIET -eq 0 ]; then
    echo ""
    echo "Successfully built:"
    for k in $BUILT_LIST; do
        echo "  - $k"
    done
fi

if [ -n "$FAILED_LIST" ]; then
    echo ""
    echo "Build failures:"
    for k in $FAILED_LIST; do
        echo "  - $k"
    done
fi

if [ -n "$SKIPPED_LIST" ] && [ $QUIET -eq 0 ]; then
    echo ""
    echo "Skipped:"
    for k in $SKIPPED_LIST; do
        echo "  - $k"
    done
fi

echo ""
echo "=============================================================================="
echo " Completed at $(date)"
echo "=============================================================================="

# Exit with error if any failed
if [ $FAILED -gt 0 ]; then
    exit 1
fi
exit 0
