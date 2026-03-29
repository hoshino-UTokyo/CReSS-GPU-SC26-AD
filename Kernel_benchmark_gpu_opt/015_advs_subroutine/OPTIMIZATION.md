# Optimization Record: advs

## Roofline Classification
- Arithmetic Intensity (DRAM): 0.26 FLOP/Byte
- Classification: Memory-bound
- Achieved BW: 1657 GB/s (41.4% of peak 4000 GB/s)
- GPU Time Rank: #7 (5.3% of total GPU time, but ~8ms per call)

## Optimizations Applied

1. **Merge independent 2nd-order flux computations**
   - Before: tmp1, tmp2, tmp3 in 3 separate `!$acc kernels` blocks.
   - After: All 3 loops in single `!$acc kernels` block (implicit sync).
   - Rationale: Eliminate 2 kernel launch overheads.

2. **Merge independent 4th-order horizontal flux computations**
   - Before: tmp1, tmp2 (4th order) in 2 separate `!$acc kernels` blocks.
   - After: Both loops in single `!$acc kernels` block.

3. **Merge vadv update + sfrc final addition (advopt==3)**
   - Before: 2 separate blocks.
   - After: Both loops in single block.

## Results
| Metric          | Before    | After     | Change  |
|-----------------|-----------|-----------|---------|
| Execution time  | 6.866 ms  | 6.832 ms  | -0.5%   |
| Speedup         | -         | 1.005x    | -       |
| Max rel error   | 0.0       | 0.0       | -       |
| Validation      | PASSED    | PASSED    | -       |

## Notes
- Minimal improvement because the advs kernel has strong inter-step data dependencies (2nd order → divergence → 4th order → update), limiting kernel fusion opportunities.
- Each step is a stencil operation reading results from the previous step, requiring synchronization.
- The 10 kernel launches (advopt=3) reduce to ~7, saving ~15us out of 6.8ms total.
- Further optimization would require algorithmic changes (e.g., operator fusion across 2nd/4th order) which changes numerical behavior.
