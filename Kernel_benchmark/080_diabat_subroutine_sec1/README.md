# Kernel 080: subroutine

## Source Location
- **File**: Src/diabat.f90
- **Subroutine**: subroutine
- **Line**: ~277
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates virtual potential temperature for dry/moist/cloud

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diabat.f90 :: subroutine s_diabat (virtual potential temperature)
! Summary : Calculates virtual potential temperature for dry/moist/cloud
!           physics cases as part of diabatic forcing computation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside this parallel region.
!   - Reads module constant epsav from comphy.
!   - No synchronization constructs.
!   - Conditional on fmois (dry/moist) and cphopt (cloud physics).
!   - All grid points are independent (embarrassingly parallel).
!   - Uses only basic arithmetic and division.
! Next:
!   - Direct OpenACC kernels should work well.
!   - Conditionals can be evaluated outside kernel for efficiency.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
