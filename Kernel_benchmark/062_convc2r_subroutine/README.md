# Kernel 062: subroutine

## Source Location
- **File**: Src/convc2r.f90
- **Subroutine**: subroutine
- **Line**: ~128

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates autoconversion rate from cloud water to rain water

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: convc2r.f90 :: subroutine s_convc2r
! Summary : Calculates autoconversion rate from cloud water to rain water
!           using threshold-based Kessler-type parameterization.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No module variable access inside loops.
!   - No synchronization constructs.
!   - Conditional logic for conversion threshold.
!   - Both qcf and qrf updated in place - no race conditions.
!   - All grid points are independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels should work well.
!   - Threshold-based conditionals may cause minor warp divergence.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
