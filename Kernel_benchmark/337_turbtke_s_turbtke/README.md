# Kernel 337: s_turbtke

## Source Location
- **File**: Src/turbtke.f90
- **Subroutine**: s_turbtke
- **Line**: ~207

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate TKE mixing term by computing divergence of turbulent

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 100.5M
- **Total Time**: 4.648s
- **Average Time per Call**: 12.912ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: turbtke.f90 :: s_turbtke
! Summary : Calculate TKE mixing term by computing divergence of turbulent
!           fluxes with terrain and map scale factor corrections
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls within loops
!   - Writes to tkefrc (in/out), tmp1, tmp2, tmp3 arrays
!   - Complex conditional structure (trnopt, mfcopt, mpopt)
!   - Multiple temporary arrays used for intermediate calculations
!   - No synchronization constructs within parallel region
! Next:
!   - GPU port requires handling multiple code paths
!   - Consider separate kernels for terrain vs non-terrain cases
!   - Temporary arrays can use shared memory or registers
!   - Map scale factor combinations may benefit from kernel specialization
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.5M
!   - TotalTime: 4.648s (0.16%)
!   - AvgTime: 12.912ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
