# Optimization Record: pgrad

## Roofline Classification
- Arithmetic Intensity (DRAM): 0.35 FLOP/Byte
- Classification: Memory-bound
- Achieved BW: 1474 GB/s (36.9% of peak 4000 GB/s)
- GPU Time Rank: #1 (11.0% of total GPU time)

## Optimizations Applied

1. **Inline tmp3 intermediate computations**
   - Before: tmp3 written 3 times as intermediate buffer (pp→tmp3, divdamp→tmp3, jcb*tmp3), then read for upg/vpg/wpg. 6 separate `!$acc kernels` blocks.
   - After: wpg, upg, vpg computed directly from pp, jcb, tmp1 without tmp3 writes. For divopt==2: `tmp3h_eff = jcb*pp + divch*tmp1` inlined into upg/vpg expressions.
   - Rationale: Eliminates ~3.6 GB of redundant DRAM traffic (3 write + 4 read passes over tmp3 array of ~394 MB each).

2. **Kernel fusion: upg + vpg in single `!$acc kernels` block**
   - Before: Separate `!$acc kernels` blocks for upg and vpg loops.
   - After: Both loops in one `!$acc kernels` block with implicit sync.
   - Rationale: Reduces kernel launch overhead and allows compiler to optimize data movement.

3. **Kernel count reduction: 6 → 2 launches (trnopt==0 path)**
   - Before: tmp3_init, wpg, tmp3_reset, tmp3_jcb, upg, vpg = 6 kernel launches
   - After: wpg (1 kernel), upg+vpg (1 kernels block with 2 loops) = 2 launches

## Results
| Metric          | Before    | After     | Change  |
|-----------------|-----------|-----------|---------|
| Execution time  | 4.026 ms  | 2.240 ms  | -44.4%  |
| Speedup         | -         | 1.80x     | -       |
| Max rel error   | 0.0       | 0.0       | -       |
| Validation      | PASSED    | PASSED    | -       |
