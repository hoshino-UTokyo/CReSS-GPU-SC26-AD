# Optimization Record: diver3d

## Roofline Classification
- Arithmetic Intensity (DRAM): 0.25 FLOP/Byte
- Classification: Memory-bound
- Achieved BW: 1682 GB/s (42.1% of peak 4000 GB/s)
- GPU Time Rank: #4 (7.7% of total GPU time)

## Optimizations Applied

1. **Inline tmp1/tmp2/tmp3 into divergence computation**
   - Before: Phase 1 writes 3 intermediate arrays (var8u*u→tmp1, rmf8v*var8v*v→tmp2, var8w*wc→tmp3), Phase 2 reads them for finite differences. 4 kernel launches.
   - After: Divergence computed directly from input arrays in single kernel. 1 kernel launch.
   - Rationale: Eliminate ~2.4 GB intermediate writes/reads.

2. **Add collapse(3)**
   - Before: Separate `!$acc loop independent` per loop dimension.
   - After: `collapse(3)` for optimal thread mapping.

## Results
| Metric          | Before    | After     | Change  |
|-----------------|-----------|-----------|---------|
| Execution time  | 2.740 ms  | 2.731 ms  | -0.3%   |
| Speedup         | -         | 1.003x    | -       |
| Max rel error   | 0.0       | 0.0       | -       |
| Validation      | PASSED    | PASSED    | -       |

## Notes
- Minimal speedup despite significant reduction in theoretical memory traffic.
- Likely cause: stencil access pattern in the inlined version reads 7 arrays at 2 stencil points each, creating high memory pressure that reduces L2 cache effectiveness.
- The original 2-phase approach benefits from temporal locality: Phase 1 writes warm L2 cache for Phase 2 reads.
- At 42% BW efficiency, further optimization would require shared memory tiling or explicit data management (beyond OpenACC kernels scope).
