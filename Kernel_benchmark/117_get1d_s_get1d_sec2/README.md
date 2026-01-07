# Kernel 117: s_get1d

## Source Location
- **File**: Src/get1d.f90
- **Subroutine**: s_get1d
- **Line**: ~849
- **Section**: 2 of 3 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Convert potential temperature to virtual potential temperature

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: get1d.f90 :: s_get1d
! Summary : Convert potential temperature to virtual potential temperature
!           using water vapor mixing ratio for moist case.
! GPU diff: Easy
! Findings:
!   - Simple 1D element-wise update on pt1d array
!   - Conditional based on gpvvar flag (outside parallel region)
!   - No reductions or synchronization
!   - Uses qv1d for moisture correction
! Next:
!   - Direct 1D GPU kernel port
!   - Small array size (4*nk-3), may be more efficient on CPU
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
