#!/bin/bash
# Binary search test runner
# Usage: bash bsearch_test.sh <disable_start_line> <disable_end_line> <run_name>
# Disables kernels from line $1 to $2 of /tmp/all_gpu_ids.txt

set -e
BASEDIR=/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload
IDS_FILE=/tmp/all_gpu_ids.txt
START=$1
END=$2
NAME=$3

if [ -z "$NAME" ]; then
  echo "Usage: $0 <start_line> <end_line> <run_name>"
  exit 1
fi

# Generate disable flags
FLAGS=$(sed -n "${START},${END}p" "$IDS_FILE" | sed 's/^/-D/' | tr '\n' ' ')
COUNT=$(sed -n "${START},${END}p" "$IDS_FILE" | wc -l)
echo "=== $NAME: disabling lines $START-$END ($COUNT kernels) ==="

# Write compile.conf
cat > "$BASEDIR/compile.conf" << EOF
################################################################################
#  $NAME: disable lines $START-$END ($COUNT kernels)
################################################################################

LDFLAGS = -fast -mp -acc -gpu=managed
FFLAGS = -fast -mp -acc -gpu=managed -Mpreprocess -Mbyteswapio -mcmodel=medium -Minfo=accel -DUSE_GPU ${FLAGS}\$(EXTRA_FLAGS)

FC      = mpifort
CC      = nvc

AR         = ar
ARFLAGS    = vru
CPP        = cpp -E
CPPFLAGS   = -DDIRECT -DIEEE -DASCII -DLEN -DUSE_GPU
AS         = as
RANLIB     = ranlib
EOF

# Clean and build
cd "$BASEDIR"
csh compile_radlib.csh clean 2>&1 | tail -1
csh compile_radlib.csh solver compile.conf 2>&1 | tail -1

# Run short simulation
cd "$BASEDIR/test_real"
rm -f result/*dmp* result/*mon* result/*geo*
mpirun -np 1 ../solver.exe < user_short.conf > "log.${NAME}.txt" 2> "log.${NAME}.err"

# Check result
DIVERGED=$(grep "ppmax" "log.${NAME}.txt" | grep -c "\-0\.9999" || true)
TOTAL_STEPS=$(grep -c "ppmax" "log.${NAME}.txt" || true)
LAST_PPMAX=$(grep "ppmax" "log.${NAME}.txt" | tail -1)

if [ "$DIVERGED" -gt 0 ]; then
  # Find divergence step
  DIV_STEP=$(grep "ppmax" "log.${NAME}.txt" | grep -n "" | grep "\-0\.9999" | head -1 | cut -d: -f1)
  echo "RESULT: **FAIL** at step ~$DIV_STEP (diverged $DIVERGED/$TOTAL_STEPS steps)"
else
  echo "RESULT: **PASS** ($TOTAL_STEPS steps completed)"
  echo "Last: $LAST_PPMAX"
fi
