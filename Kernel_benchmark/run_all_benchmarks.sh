#!/bin/bash
#==============================================================================
# Run All Kernel Benchmarks
#==============================================================================
# Usage: ./run_all_benchmarks.sh [OPTIONS]
#
# Options:
#   -j N        Number of OpenMP threads (default: 4)
#   -q          Quiet mode (only show summary)
#   -s          Stop on first failure
#   -h          Show this help
#==============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Default settings
OMP_THREADS=72
QUIET=0
STOP_ON_FAIL=0

# Parse arguments
while getopts "j:qsh" opt; do
    case $opt in
        j) OMP_THREADS=$OPTARG ;;
        q) QUIET=1 ;;
        s) STOP_ON_FAIL=1 ;;
        h)
            head -15 "$0" | tail -14
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
PASSED=0
FAILED=0
SKIPPED=0
NO_BENCHMARK=0

# Results storage
PASSED_LIST=""
FAILED_LIST=""

# Header
echo "=============================================================================="
echo " Kernel Benchmark Runner"
echo "=============================================================================="
echo " Directory: $SCRIPT_DIR"
echo " OpenMP Threads: $OMP_THREADS"
echo " Date: $(date)"
echo "=============================================================================="
echo ""

export OMP_NUM_THREADS=$OMP_THREADS
export OMP_STACKSIZE=256000

# Iterate through all kernel directories
for dir in [0-9][0-9][0-9]_*/; do
    [ -d "$dir" ] || continue

    TOTAL=$((TOTAL + 1))
    kernel_name=$(basename "$dir")

    # Check if benchmark executable exists
    if [ ! -f "$dir/kernel_benchmark" ]; then
        # Check if source exists but not compiled
        if [ -f "$dir/kernel_benchmark.f90" ]; then
            if [ $QUIET -eq 0 ]; then
                echo "[$kernel_name] SKIP: Not compiled (run 'make' in directory)"
            fi
            SKIPPED=$((SKIPPED + 1))
        else
            if [ $QUIET -eq 0 ]; then
                echo "[$kernel_name] SKIP: No benchmark available"
            fi
            NO_BENCHMARK=$((NO_BENCHMARK + 1))
        fi
        continue
    fi

    # Check if data directory exists
    if [ ! -d "$dir/data" ]; then
        if [ $QUIET -eq 0 ]; then
            echo "[$kernel_name] SKIP: No data directory"
        fi
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Check if data directory has actual data files (not just params.txt)
    # Use -L to follow symlinks
    bin_files=$(find -L "$dir/data" -maxdepth 1 -name "*.bin" 2>/dev/null | head -1)
    if [ -z "$bin_files" ]; then
        if [ $QUIET -eq 0 ]; then
            echo "[$kernel_name] SKIP: No data files (dump incomplete)"
        fi
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Run benchmark
    if [ $QUIET -eq 0 ]; then
        echo ""
        echo "----------------------------------------------------------------------"
        echo "[$kernel_name] Running..."
        echo "----------------------------------------------------------------------"
    fi

    cd "$dir"

    # Execute benchmark and capture output
    output=$(./kernel_benchmark 2>&1)
    exit_code=$?

    cd "$SCRIPT_DIR"

    # Check result
    if [ $exit_code -eq 0 ]; then
        PASSED=$((PASSED + 1))
        PASSED_LIST="$PASSED_LIST $kernel_name"
        if [ $QUIET -eq 0 ]; then
            echo "$output" | grep -E "(Average time|Validation)"
            echo "[$kernel_name] PASSED"
        fi
    else
        FAILED=$((FAILED + 1))
        FAILED_LIST="$FAILED_LIST $kernel_name"
        if [ $QUIET -eq 0 ]; then
            echo "$output" | grep -E "(Average time|Max relative error|Error count|Validation)"
            echo "[$kernel_name] FAILED"
        fi

        if [ $STOP_ON_FAIL -eq 1 ]; then
            echo ""
            echo "Stopping on first failure."
            break
        fi
    fi
done

# Summary
echo ""
echo "=============================================================================="
echo " Summary"
echo "=============================================================================="
echo " Total directories:    $TOTAL"
echo " Passed:               $PASSED"
echo " Failed:               $FAILED"
echo " Skipped (not built):  $SKIPPED"
echo " No benchmark:         $NO_BENCHMARK"
echo "=============================================================================="

if [ -n "$PASSED_LIST" ]; then
    echo ""
    echo "Passed benchmarks:"
    for k in $PASSED_LIST; do
        echo "  - $k"
    done
fi

if [ -n "$FAILED_LIST" ]; then
    echo ""
    echo "Failed benchmarks:"
    for k in $FAILED_LIST; do
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
