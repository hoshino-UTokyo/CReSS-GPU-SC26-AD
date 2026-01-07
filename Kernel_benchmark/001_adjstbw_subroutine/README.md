# Kernel 001: subroutine

## Source Location
- **File**: Src/adjstbw.f90
- **Subroutine**: subroutine
- **Line**: ~117

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Adjusts mean water mass for bin microphysics to stay within

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: adjstbw.f90 :: subroutine s_adjstbw
! Summary : Adjusts mean water mass for bin microphysics to stay within
!           bin boundaries, resetting concentration if mass is invalid.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Pure arithmetic and conditional logic only.
!   - All grid points independent (embarrassingly parallel).
!   - 4D array access with bin categories in outer loop.
!   - Uses intrinsic conditionals, no sync constructs.
! Next:
!   - Direct OpenACC kernels with collapse for (n,k,j,i).
!   - bmw array is small and can be copied to device.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
