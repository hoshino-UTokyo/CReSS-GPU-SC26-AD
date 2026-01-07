# Kernel 312: s_steptund

## Source Location
- **File**: Src/steptund.f90
- **Subroutine**: s_steptund
- **Line**: ~577
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Set bottom boundary condition and convert Celsius to Kelvin

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: steptund.f90 :: s_steptund
! Summary : Set bottom boundary condition and convert Celsius to Kelvin
!           for soil/sea temperature output
! GPU diff: Easy
! Findings:
!   - Simple conditional for sfcopt and land type
!   - Straightforward array update with constant offset
!   - No function calls within parallel region
! Next:
!   - Collapse loops for GPU parallelization
!   - Keep tundf array resident on GPU from previous kernel
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
