# Kernel 015: subroutine

## Source Location
- **File**: Src/advs.f90
- **Subroutine**: subroutine
- **Line**: ~207

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates scalar variable advection using 2nd or 4th order

## Runtime Profile (from test_real)
- **Calls**: 3960
- **Average Loop Length**: 100.5M
- **Total Time**: 155.185s
- **Average Time per Call**: 39.188ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: advs.f90 :: subroutine s_advs
! Summary : Calculates scalar variable advection using 2nd or 4th order
!           centered finite difference schemes with mass-weighted fluxes.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constants (oned24, fourd3) from commath.
!   - No synchronization constructs.
!   - Multiple code paths based on advopt (1=2nd order full, 2=2nd+4th,
!     3=2nd+4th+vadv separated).
!   - Multi-stage stencil with temporary arrays (tmp1, tmp2, tmp3, vadv).
!   - All loops are embarrassingly parallel within each stage.
! Next:
!   - Split into kernels matching the loop structure.
!   - Temporary arrays already allocated - good for GPU data management.
!   - Consider fusing stages for reduced memory traffic.
! Runtime:
!   - Calls: 3960
!   - AvgLoops: 100.5M
!   - TotalTime: 155.185s (5.21%)
!   - AvgTime: 39.188ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
