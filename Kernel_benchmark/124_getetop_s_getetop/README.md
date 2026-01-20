# Kernel 124: s_getetop

## Source Location
- **File**: Src/getetop.f90
- **Subroutine**: s_getetop
- **Line**: ~165

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Computes z-coordinates at scalar points, radar echo top height,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getetop.f90 :: s_getetop
! Summary : Computes z-coordinates at scalar points, radar echo top height,
!           and total precipitation mixing ratio from radar hydrometeor data
! GPU diff: Medium
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic max() function - GPU compatible
!   - Multiple conditional branches based on ngropt and haiopt
!   - Potential race condition on etop(i,j) with max() update across k-loop
!   - Module variables lim34n, lim35n, lim36n, qpmin used from m_commath
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Handle etop update carefully - may need atomic or reduction approach
!   - Consider collapsing i,j loops for better GPU parallelism
!   - Ensure module constants are accessible on device
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
