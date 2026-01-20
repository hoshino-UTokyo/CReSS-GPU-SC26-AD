# Kernel 205: s_move2d

## Source Location
- **File**: Src/move2d.f90
- **Subroutine**: s_move2d
- **Line**: ~116

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Subtract grid moving velocity (umove, vmove) from horizontally

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: move2d.f90 :: s_move2d
! Summary : Subtract grid moving velocity (umove, vmove) from horizontally
!           averaged velocity profiles u1d and v1d
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple subtraction operation on 1D arrays
!   - Writes to u1d, v1d arrays (in-place modification)
!   - Single worksharing construct with 1D loop
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Small array size (nlev) may not benefit significantly from GPU offloading
!   - Consider keeping on CPU if nlev is small
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
