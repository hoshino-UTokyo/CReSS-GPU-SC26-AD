# Kernel 008: subroutine

## Source Location
- **File**: Src/adjstqa.f90
- **Subroutine**: subroutine
- **Line**: ~105

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Forces aerosol mixing ratio to be non-negative using max(0).

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: adjstqa.f90 :: subroutine s_adjstqa
! Summary : Forces aerosol mixing ratio to be non-negative using max(0).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Pure max() operation, GPU compatible.
!   - 4D array with aerosol categories in outer loop.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse for (n,k,j,i).
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
