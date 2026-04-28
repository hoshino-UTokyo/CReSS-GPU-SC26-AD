#!/bin/bash
# Quick check script for binary search tests
# Usage: ./check_nan.sh <logfile>
# Checks t=265s and t=270s for NaN (sentinel value -0.99999996E+36)

if [ -z "$1" ]; then
    echo "Usage: $0 <logfile>"
    echo "Available logs:"
    ls -lt test_real/log.solver.*.txt 2>/dev/null | head -10
    exit 1
fi

LOG="$1"

if [ ! -f "$LOG" ]; then
    echo "Error: $LOG not found"
    exit 1
fi

echo "=== Checking $LOG ==="
echo ""

# Check last timestep
LAST_TIME=$(grep "time," "$LOG" | tail -1)
echo "Last timestep: $LAST_TIME"
echo ""

# Check t=265s (should be normal in all cases)
echo "--- t=265s ---"
grep -A 8 "time, 265.000" "$LOG" | grep -E "umax|vmax|wmax|ppmax"
echo ""

# Check t=270s (NaN occurs here in GPU version)
echo "--- t=270s ---"
grep -A 8 "time, 270.000" "$LOG" | grep -E "umax|vmax|wmax|ppmax"
echo ""

# Check for sentinel values
if grep -q "\-0.99999996E+36" "$LOG"; then
    echo "*** NaN DETECTED (sentinel -0.99999996E+36 found) ***"
    echo "First occurrence:"
    grep -m 1 -B 5 "\-0.99999996E+36" "$LOG"
else
    echo "*** NO NaN detected - test PASSED ***"
fi
