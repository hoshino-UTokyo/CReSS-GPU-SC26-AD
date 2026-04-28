#!/bin/bash
# Single kernel enable test
# Usage: bash bsearch_single.sh <kernel_id> <run_name>
# Starts from baseline (all 118-127 disabled), enables ONE kernel

set -e
BASEDIR=/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload
ENABLE_ID=$1
NAME=$2

if [ -z "$NAME" ]; then
  echo "Usage: $0 <kernel_id_to_enable> <run_name>"
  exit 1
fi

# All 10 suspect kernels
ALL_IDS="294 295 296 302 303 304 305 308 310 311"

# Build disable flags, skipping the one we want to enable
FLAGS=""
for id in $ALL_IDS; do
  if [ "$id" != "$ENABLE_ID" ]; then
    FLAGS="$FLAGS -DDISABLE_GPU_${id}"
  fi
done

echo "=== $NAME: enabling only GPU_$ENABLE_ID (9 others disabled) ==="

# Write compile.conf
cat > "$BASEDIR/compile.conf" << EOF
################################################################################
#  $NAME: enable GPU_$ENABLE_ID only (9 others disabled)
################################################################################

LDFLAGS = -fast -mp -acc -gpu=managed
FFLAGS = -fast -mp -acc -gpu=managed -Mpreprocess -Mbyteswapio -mcmodel=medium -Minfo=accel -DUSE_GPU${FLAGS} \$(EXTRA_FLAGS)

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
  DIV_STEP=$(grep "ppmax" "log.${NAME}.txt" | grep -n "" | grep "\-0\.9999" | head -1 | cut -d: -f1)
  echo "RESULT: **FAIL** at step ~$DIV_STEP (diverged $DIVERGED/$TOTAL_STEPS steps)"
else
  echo "RESULT: **PASS** ($TOTAL_STEPS steps completed)"
  echo "Last: $LAST_PPMAX"
fi
