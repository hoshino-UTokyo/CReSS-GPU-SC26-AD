# Kernel 118: s_get1d

## Source Location
- **File**: Src/get1d.f90
- **Subroutine**: s_get1d
- **Line**: ~937
- **Section**: 3 of 3 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Recalculate potential temperature from virtual and convert

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: get1d.f90 :: s_get1d
! Summary : Recalculate potential temperature from virtual and convert
!           Exner function to pressure using exponential.
! GPU diff: Easy
! Findings:
!   - Two 1D loops: pt1d recalculation and p1d conversion
!   - Uses exp and log intrinsic functions
!   - No reductions or synchronization
!   - Small array size (4*nk-3)
! Next:
!   - Direct 1D GPU kernel port
!   - Consider keeping on CPU due to small problem size
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
