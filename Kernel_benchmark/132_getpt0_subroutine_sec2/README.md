# Kernel 132: subroutine

## Source Location
- **File**: Src/getpt0.f90
- **Subroutine**: subroutine
- **Line**: ~483
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets sine-curved initial potential temperature perturbation

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getpt0.f90 :: subroutine s_getpt0 (sine perturbation)
! Summary : Sets sine-curved initial potential temperature perturbation
!           for wave-like thermal initialization experiments.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constant cc from commath.
!   - No synchronization constructs.
!   - Uses intrinsic sin(), cos() - GPU compatible.
!   - Height-based conditional for perturbation region.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i) loops.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
