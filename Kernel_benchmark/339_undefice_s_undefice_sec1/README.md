# Kernel 339: s_undefice

## Source Location
- **File**: Src/undefice.f90
- **Subroutine**: s_undefice
- **Line**: ~150
- **Section**: 1 of 3 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Count valid sea ice data points within valid range (0-100)

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: undefice.f90 :: s_undefice
! Summary : Count valid sea ice data points within valid range (0-100)
!           using parallel reduction for error checking
! GPU diff: Easy
! Findings:
!   - Uses reduction(+: rstat) for counting valid points
!   - Simple 2D loop with no function calls
!   - No global writes, only local reduction
! Next:
!   - Convert to OpenACC with reduction clause
!   - Data should be present on GPU from caller
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
