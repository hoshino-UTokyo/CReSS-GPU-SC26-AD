#!/bin/bash
# TKE binary search: disable a range of suspect kernels to find TKE bug
# Usage: bash bsearch_tke.sh <disable_ids_comma_sep> <run_name>
# Example: bash bsearch_tke.sh "312,313,315,318,320" tke_half1

set -e
BASEDIR=/work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload
DISABLE_IDS=$1
NAME=$2

if [ -z "$NAME" ]; then
  echo "Usage: $0 <disable_ids_comma_sep> <run_name>"
  exit 1
fi

# Build disable flags
FLAGS=""
IFS=',' read -ra IDS <<< "$DISABLE_IDS"
for id in "${IDS[@]}"; do
  FLAGS="$FLAGS -DDISABLE_GPU_${id}"
done

echo "=== $NAME: disabling GPU_${DISABLE_IDS} ==="

# Write compile.conf
cat > "$BASEDIR/compile.conf" << EOF
################################################################################
#  $NAME: disable GPU_${DISABLE_IDS}
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

# Check TKE at step 2
TKE_STEP2=$(grep "tkemax" "log.${NAME}.txt" | sed -n '2p' | awk '{print $3}')
TOTAL_STEPS=$(grep -c "ppmax" "log.${NAME}.txt" || true)
DIVERGED=$(grep "ppmax" "log.${NAME}.txt" | grep -c "\-0\.9999" || true)

echo "TKE at step 2: $TKE_STEP2"
if [ "$DIVERGED" -gt 0 ]; then
  echo "RESULT: DIVERGED at step ~$(grep "ppmax" "log.${NAME}.txt" | grep -n "" | grep "\-0\.9999" | head -1 | cut -d: -f1)"
elif [ "$TKE_STEP2" = "0.00000000E+00" ]; then
  echo "RESULT: TKE_PASS (step2 TKE=0, correct)"
else
  echo "RESULT: TKE_FAIL (step2 TKE non-zero)"
fi
