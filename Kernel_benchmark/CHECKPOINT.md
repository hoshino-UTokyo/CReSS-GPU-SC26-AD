# Kernel Benchmark Progress Checkpoint

**Date**: 2026-01-07
**Status**: Complete - All 20 Benchmarks Passing

## Latest Update (2026-01-07)

### Work Completed
1. **Updated benchmark_instruction.md** with comprehensive data dump rules:
   - Added "CRITICAL RULE: Dump ALL variables used inside the OpenMP parallel section"
   - Documented categories: Arrays, Scalar Parameters, Derived Constants, Module Constants
   - Added "Common Mistakes to Avoid" section

2. **Fixed 230_pgrad benchmark**:
   - Instead of regenerating dump data, modified benchmark to compute tmp1 using diver3d
   - Added diver3d kernel to kernel_benchmark.f90
   - Reads all diver3d inputs (mf, rmf, rmf8u, rmf8v, rst8u, rst8v, rst8w, u, v, wc)
   - Computes divergence (tmp1) before running pgrad kernel
   - Result: Max error = 0, perfect match with reference

3. **All 20 benchmarks now passing**:
   - 012_advbspi_subroutine, 015_advs_subroutine, 027_bc4news_s_bc4news
   - 038_bcycle_s_bcycle, 067_copy3d_s_copy3d, 083_diagni_s_diagni
   - 087_disptke_s_disptke, 090_diver2d_s_diver2d, 091_diver3d_s_diver3d
   - 094_diverpiv_s_diverpiv, 230_pgrad_subroutine, 236_phy2cnt_s_phy2cnt
   - 308_steppi_s_steppi, 310_steps_s_steps, 313_stepuv_subroutine
   - 315_stepwi_subroutine, 361_vbcu_s_vbcu, 362_vbcv_s_vbcv
   - 363_vbcw_s_vbcw, 364_vbcwc_s_vbcwc

### Key Fix: pgrad tmp1 Issue
The pgrad kernel requires tmp1 (divergence) when divopt >= 1. Rather than regenerating dump data:
- Added kernel_diver3d to the benchmark
- Reads necessary arrays from dump data (already present in pgrad dump)
- Computes tmp1 at benchmark runtime
- This approach is more robust: works without needing to regenerate dump data

### Files Modified
- `benchmark_instruction.md`: Added dump rules
- `Src/pgrad.f90`: Added tmp1 dump (line 353-356)
- `Kernel_benchmark/230_pgrad_subroutine/kernel_benchmark.f90`: Added tmp1_input handling

---

## Previous Status

## Active Job
- **Job ID**: 1256130
- **Status**: QUEUED
- **Queue**: short-g
- **Estimated Time**: 2.5 hours
- **Purpose**: Generate dump data for remaining kernels

## Current Progress

### Phase 1: Data Dump
- **Completed**: 20 kernels have dump data in `test_real/kernel_dump/`
- **Dump code added**: 158 source files have dump instrumentation
- **Need simulation run**: ~130 kernels still need dump data generated

### Phase 2: Benchmark Creation
- **Completed**: 30 kernel benchmarks
- **All kernels with dump data have benchmarks created**

## Completed Benchmarks (30)
1. 012_advbspi_subroutine
2. 014_advp_subroutine
3. 015_advs_subroutine
4. 016_advuvw_subroutine
5. 027_bc4news_s_bc4news
6. 038_bcycle_s_bcycle
7. 041_bruntv_s_bruntv
8. 046_buoywsi_s_buoywsi
9. 067_copy3d_s_copy3d
10. 073_curveuvw_s_curveuvw
11. 083_diagni_s_diagni
12. 087_disptke_s_disptke
13. 090_diver2d_s_diver2d
14. 091_diver3d_s_diver3d
15. 094_diverpiv_s_diverpiv
16. 115_gaussel_subroutine
17. 230_pgrad_subroutine
18. 231_pgradiv_s_pgradiv
19. 236_phy2cnt_s_phy2cnt
20. 303_smoo4s_subroutine
21. 304_smoo4uvw_s_smoo4uvw
22. 308_steppi_s_steppi
23. 310_steps_s_steps
24. 313_stepuv_subroutine
25. 315_stepwi_subroutine
26. 336_turbs_subroutine
27. 361_vbcu_s_vbcu
28. 362_vbcv_s_vbcv
29. 363_vbcw_s_vbcw
30. 364_vbcwc_s_vbcwc

## Kernels with Dump Data (20)
Located in `test_real/kernel_dump/`:
- advbspi, advs, bc4news, bcycle, buoywsi
- copy3d, diagni, disptke, diver2d, diver3d
- pgrad, phy2cnt, steppi, steps, stepuv
- stepwi, vbcu, vbcv, vbcw, vbcwc

## Resume Instructions

When resuming work:

1. **Check job status**:
   ```bash
   qstat 1256130
   ```

2. **If job completed**, check for new dump data:
   ```bash
   ls test_real/kernel_dump/
   ```

3. **Copy new dump data** to Kernel_benchmark directories:
   ```bash
   for k in $(ls test_real/kernel_dump/); do
     dest=$(ls -d Kernel_benchmark/*_${k}_* 2>/dev/null | head -1)
     if [ -n "$dest" ] && [ ! -d "$dest/data" ]; then
       mkdir -p "$dest/data"
       cp -r test_real/kernel_dump/$k/* "$dest/data/"
       echo "Copied: $k -> $dest"
     fi
   done
   ```

4. **Create benchmarks** for new kernels following the template pattern.

## Next Steps

### Step 1: Generate More Dump Data
Run simulation to generate dump data for remaining kernels:
```bash
cd /work/jh250015/g24000/SC26/CReSS3.5.1m_SPN_RAD1.4.3_20230323_upload/test_real
qsub MIYABI-solver.sh
```

Wait for job to complete and check `kernel_dump/` for new data.

### Step 2: Copy Dump Data
For each new kernel dump, copy to Kernel_benchmark:
```bash
cp -r kernel_dump/<kernel_name>/* ../Kernel_benchmark/<id>_<kernel_name>_<subroutine>/data/
```

### Step 3: Create Benchmarks
For each kernel with new dump data:
1. Read source file and README
2. Create kernel_benchmark.f90 following template from 012_advbspi_subroutine
3. Create Makefile
4. Test build and run

## Priority Kernels (High Runtime)

Based on omp_profile.txt, these kernels have high execution time:
1. stepwi (410.58s) - DONE
2. pgrad (296.86s) - DONE
3. diver3d (232.45s) - DONE
4. advbspi (180.04s) - DONE
5. diver2d (158.99s) - DONE
6. stepuv (158.55s) - DONE
7. advs (151.49s) - DONE
8. buoywsi (143.81s) - DONE
9. gaussel (141.52s) - DONE
10. pgradiv (115.49s) - DONE

Next priority (not yet done):
11. diverpiv (105.20s) - DONE
12. smoo4s (88.88s) - DONE
13. smoo4uvw (28.61s) - DONE
14. newblk (16.82s) - TODO: needs dump data
15. curveuvw (17.31s) - DONE
16. upwqp (17.46s) - TODO: needs dump data

## Technical Notes

- Dump module: `Kernel_benchmark/dump_kernel_data.f90` (compiled)
- Solver executable includes dump functionality
- Each kernel dumps at its final call (DUMP_TARGET_CALL from omp_profile.txt Count)
- Binary format: stream access, unformatted, column-major order

## Files Created for This Session

- `Kernel_benchmark/find_incomplete.sh` - Script to find incomplete kernels
- `Kernel_benchmark/check_data.sh` - Script to check data copy status
- `Kernel_benchmark/CHECKPOINT.md` - This checkpoint file
