# Kernel 383: subroutine

## Source Location
- **File**: Src/vsps.f90
- **Subroutine**: subroutine
- **Line**: ~162

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Applies vertical sponge damping near model top for scalar

## Runtime Profile (from test_real)
- **Calls**: 2160
- **Average Loop Length**: 101.2M
- **Total Time**: 12.669s
- **Average Time per Call**: 5.865ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vsps.f90 :: subroutine s_vsps
! Summary : Applies vertical sponge damping near model top for scalar
!           variables, relaxing toward GPV data or base state.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Simple arithmetic with rbct damping coefficients.
!   - Only upper levels computed (k >= ksp0 - sparse in k).
!   - Conditional on vspopt for GPV vs base state damping target.
! Next:
!   - Direct OpenACC kernels should work well.
!   - Upper-level-only computation - consider k-range optimization.
! Runtime:
!   - Calls: 2160
!   - AvgLoops: 101.2M
!   - TotalTime: 12.669s (0.43%)
!   - AvgTime: 5.865ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
