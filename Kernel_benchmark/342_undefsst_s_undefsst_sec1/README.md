# Kernel 342: s_undefsst

## Source Location
- **File**: Src/undefsst.f90
- **Subroutine**: s_undefsst
- **Line**: ~155
- **Section**: 1 of 3 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Count valid SST data points within valid range (268.16-323.16 K)

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: undefsst.f90 :: s_undefsst
! Summary : Count valid SST data points within valid range (268.16-323.16 K)
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
