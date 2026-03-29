# Optimization Record: gaussel

## Roofline Classification
- Arithmetic Intensity (DRAM): 0.44 FLOP/Byte
- Classification: Memory-bound / Latency-bound (many small kernel launches)
- Achieved BW: 812-1503 GB/s (20-38% of peak)
- GPU Time Rank: #3 (8.1% of total GPU time)
- CPU→GPU speedup: 7.7x (lowest among top 10)

## Optimizations Applied

1. **Consolidate k-sequential loops into single `!$acc kernels` blocks**
   - Before: Each k-level in the Thomas algorithm (forward elimination + back substitution) was in a separate `!$acc kernels` block with `!$acc wait`. With kend-kstr ≈ 123 levels, this produced ~247 kernel launches, each with ~5-10us host-side overhead.
   - After: Forward elimination (first level + middle loop + final level) in 1 `!$acc kernels` block. Back substitution in 1 `!$acc kernels` block. Total: 2 blocks instead of ~247.
   - Rationale: Kernel launch overhead was ~40% of total execution time (~1.9ms out of 4.67ms).

2. **Add collapse(2) to j,i loops**
   - Before: Separate `!$acc loop independent` per dimension.
   - After: `!$acc loop independent collapse(2)` for optimal 2D thread mapping.

3. **Remove redundant `!$acc wait` calls**
   - Before: `!$acc wait` after each per-k kernel launch.
   - After: Implicit synchronization within `!$acc kernels` block handles ordering.

## Results
| Metric          | Before    | After     | Change  |
|-----------------|-----------|-----------|---------|
| Execution time  | 4.674 ms  | 2.556 ms  | -45.3%  |
| Speedup         | -         | 1.83x     | -       |
| Max rel error   | 0.0       | 0.0       | -       |
| Validation      | PASSED    | PASSED    | -       |

## Notes
- This optimization reduces HOST-SIDE overhead, not GPU kernel efficiency.
- The underlying Thomas algorithm is still k-sequential on the GPU.
- Further improvement would require Parallel Cyclic Reduction (PCR) or similar parallel tridiagonal solver, which is a major algorithmic change.
