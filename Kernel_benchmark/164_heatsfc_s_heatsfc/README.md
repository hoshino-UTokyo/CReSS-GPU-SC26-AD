# Kernel 164: s_heatsfc

## Source Location
- **File**: Src/heatsfc.f90
- **Subroutine**: s_heatsfc
- **Line**: ~158

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Compute sensible and latent heat fluxes on surface based on

## Runtime Profile (from test_real)
- **Calls**: 361
- **Average Loop Length**: 806.4K
- **Total Time**: 0.058s
- **Average Time per Call**: 0.161ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: heatsfc.f90 :: s_heatsfc
! Summary : Compute sensible and latent heat fluxes on surface based on
!           land type, temperature, and moisture conditions.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls; uses intrinsic exp/log
!   - Complex branching based on fmois (dry/moist) and land type
!   - Reads from t, qv, qvsfc, ct, cq, kai, tund, tice
!   - Writes to hs, le arrays (output)
!   - No sync constructs
! Next:
!   - Convert to OpenACC with collapse(2) on j,i loops
!   - Branch logic based on land type may cause GPU thread divergence
!   - Consider separating dry/moist cases into different kernels
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.058s (0.00%)
!   - AvgTime: 0.161ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
