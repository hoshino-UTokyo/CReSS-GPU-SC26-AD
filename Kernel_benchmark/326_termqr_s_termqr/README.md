# Kernel 326: s_termqr

## Source Location
- **File**: Src/termqr.f90
- **Subroutine**: s_termqr
- **Line**: ~132

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate terminal velocity of rain water using power-law

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: termqr.f90 :: s_termqr
! Summary : Calculate terminal velocity of rain water using power-law
!           formulation based on density and mixing ratio
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Calls intrinsic functions only (exp, log, sqrt)
!   - Single output array (urq)
!   - Simple conditional (threshold check)
!   - No synchronization constructs
! Next:
!   - Straightforward GPU port with OpenACC or OpenACC
!   - Single kernel (data managed via Unified Memory)
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
