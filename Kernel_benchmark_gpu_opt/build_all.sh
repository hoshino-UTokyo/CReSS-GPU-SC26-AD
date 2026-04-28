#!/bin/bash
#
# Build all GPU kernel benchmarks
#

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "========================================"
echo " Building all GPU Kernel Benchmarks"
echo "========================================"

SUCCESS=0
FAILED=0
FAILED_LIST=""

for dir in */; do
    if [ -f "${dir}Makefile" ]; then
        echo ""
        echo "--- Building: ${dir%/} ---"
        cd "$dir"
        if make clean > /dev/null 2>&1 && make; then
            echo "[OK] ${dir%/}"
            SUCCESS=$((SUCCESS + 1))
        else
            echo "[FAILED] ${dir%/}"
            FAILED=$((FAILED + 1))
            FAILED_LIST="$FAILED_LIST ${dir%/}"
        fi
        cd "$SCRIPT_DIR"
    fi
done

echo ""
echo "========================================"
echo " Build Summary"
echo "========================================"
echo " Success: $SUCCESS"
echo " Failed:  $FAILED"
if [ -n "$FAILED_LIST" ]; then
    echo " Failed benchmarks:$FAILED_LIST"
fi
echo "========================================"

exit $FAILED
