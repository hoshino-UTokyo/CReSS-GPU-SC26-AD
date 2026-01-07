# Kernel 262: s_rotuvm2s

## Source Location
- **File**: Src/rotuvm2s.f90
- **Subroutine**: s_rotuvm2s
- **Line**: ~194

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Rotate velocity components from projected grid to lat/lon grid for different map projections

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rotuvm2s.f90 :: s_rotuvm2s
! Summary : Rotate velocity components from projected grid to lat/lon grid for different map projections
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic functions: cos, sin
!   - Three conditional branches based on mpopt (1, 2, or 4)
!   - Outer k loop is serial with inner !$omp do on i,j
!   - In-place update of u and v arrays using temp variables utmp, vtmp
!   - Conditional on u,v > lim34n before rotation
!   - No synchronization constructs
! Next:
!   - Collapse k,j,i loops for better GPU parallelism
!   - Consider separate kernels for each mpopt case
!   - Use OpenACC teams distribute parallel for collapse(3)
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
