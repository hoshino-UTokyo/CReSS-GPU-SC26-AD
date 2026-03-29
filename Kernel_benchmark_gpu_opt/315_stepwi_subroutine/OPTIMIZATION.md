# Optimization Record: stepwi

## Roofline Classification
- Arithmetic Intensity (DRAM): 0.60 FLOP/Byte
- Classification: Memory-bound
- Achieved BW: 2151 GB/s (53.8% of peak 4000 GB/s)
- GPU Time Rank: #2 (11.0% of total GPU time)

## Optimizations Applied

1. **Inline tmp1/tmp2 intermediates in Step 2**
   - Before: Step 2a writes tmp1=g05*rbr/rcsq and tmp2=dziv/jcb to arrays. Step 2b reads them to compute fw=tmp1+tmp2, wc=tmp1-tmp2. 2 separate kernels.
   - After: Scalar temporaries val1/val2 replace array writes. fw and wc computed in single loop from scalars.
   - Rationale: Eliminates 2×~420MB intermediate array writes + reads.

2. **Add collapse(3) to all loops**
   - Before: Separate `!$acc loop independent` per dimension (3 directives per loop nest).
   - After: Single `!$acc loop independent collapse(3)` per loop nest.

3. **Kernel fusion: Step 1 + Step 2 in single `!$acc kernels` block**
   - Before: 5 separate kernel launches.
   - After: 2 kernel launches (Step 1+2 fused, Step 3 separate).

## Results
| Metric          | Before    | After     | Change  |
|-----------------|-----------|-----------|---------|
| Execution time  | 3.885 ms  | 3.210 ms  | -17.4%  |
| Speedup         | -         | 1.21x     | -       |
| Max rel error   | 1.43e-7   | 1.43e-7   | -       |
| Validation      | PASSED    | PASSED    | -       |
