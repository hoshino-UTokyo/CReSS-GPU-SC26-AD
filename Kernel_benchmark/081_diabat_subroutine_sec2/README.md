# Kernel 081: subroutine

## Source Location
- **File**: Src/diabat.f90
- **Subroutine**: subroutine
- **Line**: ~429
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Adds time tendency to diabatic term and computes final diabatic

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diabat.f90 :: subroutine s_diabat (diabatic value computation)
! Summary : Adds time tendency to diabatic term and computes final diabatic
!           value for pressure equation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Two stages: (1) add time tendency, (2) normalize by ptv.
!   - All grid points are independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels for each loop nest.
!   - Consider fusing the two stages into single kernel.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
