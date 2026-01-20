# Kernel 238: s_phycood

## Source Location
- **File**: Src/phycood.f90
- **Subroutine**: s_phycood
- **Line**: ~290
- **Section**: 2 of 3 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Find the lowest flat level index using min reduction over k levels.

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: phycood.f90 :: s_phycood
! Summary : Find the lowest flat level index using min reduction over k levels.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses reduction(min: kflat) clause
!   - Reads from zsth array, no writes to shared arrays
!   - Simple 1D loop with conditional and reduction
! Next:
!   - Use OpenACC parallel loop with reduction(min:kflat)
!   - Small loop range (nk typically ~50-100), may be better on CPU
!   - Consider keeping this on host if nk is small
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
