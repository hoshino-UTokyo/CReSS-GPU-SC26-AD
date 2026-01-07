# Kernel 263: s_rotuvs2m

## Source Location
- **File**: Src/rotuvs2m.f90
- **Subroutine**: s_rotuvs2m
- **Line**: ~169

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Rotate velocity components from lat/lon grid to projected model grid for different projections

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rotuvs2m.f90 :: s_rotuvs2m
! Summary : Rotate velocity components from lat/lon grid to projected model grid for different projections
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic functions: cos, sin, tan
!   - Three conditional branches based on mpopt (1, 2, or 4)
!   - Outer kd loop is serial with inner !$omp do on id,jd
!   - In-place update of udat and vdat arrays using temp variables
!   - Conditional on udat,vdat > lim34n before rotation
!   - No synchronization constructs
! Next:
!   - Collapse kd,jd,id loops for better GPU parallelism
!   - Consider separate kernels for each mpopt case
!   - Use OpenACC teams distribute parallel for collapse(3)
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
