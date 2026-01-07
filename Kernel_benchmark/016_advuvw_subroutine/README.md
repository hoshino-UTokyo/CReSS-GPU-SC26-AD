# Kernel 016: subroutine

## Source Location
- **File**: Src/advuvw.f90
- **Subroutine**: subroutine
- **Line**: ~231

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates velocity (u, v, w) advection using 2nd or 4th order

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 100.6M
- **Total Time**: 43.681s
- **Average Time per Call**: 121.336ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: advuvw.f90 :: subroutine s_advuvw
! Summary : Calculates velocity (u, v, w) advection using 2nd or 4th order
!           centered finite difference schemes with mass-weighted fluxes.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constants (oned24, fourd3) from commath.
!   - No synchronization constructs.
!   - Multiple code paths based on advopt (1=2nd order, 2=2nd+4th, 3=separated).
!   - Multi-stage stencil with temporary arrays (tmp1, tmp2, tmp3, hadv, vadv).
!   - Three separate velocity components (u, v, w) computed sequentially.
!   - All loops are embarrassingly parallel within each stage.
! Next:
!   - Can use multiple OpenACC kernels matching the loop structure.
!   - Consider computing u, v, w advection in parallel if temporary arrays
!     are independent.
!   - Temporary arrays already allocated - good for GPU data management.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.6M
!   - TotalTime: 43.681s (1.47%)
!   - AvgTime: 121.336ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
