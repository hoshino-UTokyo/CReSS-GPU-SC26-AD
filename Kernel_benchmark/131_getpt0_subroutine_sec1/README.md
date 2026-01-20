# Kernel 131: subroutine

## Source Location
- **File**: Src/getpt0.f90
- **Subroutine**: subroutine
- **Line**: ~319
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets bubble-shaped initial potential temperature perturbation

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getpt0.f90 :: subroutine s_getpt0 (bubble perturbation)
! Summary : Sets bubble-shaped initial potential temperature perturbation
!           using cosine-squared profile for thermal bubble experiments.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constant cc from commath.
!   - No synchronization constructs.
!   - Uses intrinsic sqrt(), cos() - GPU compatible.
!   - Outer loop over bubble number (pt0num, typically small).
!   - Conditional update based on distance from bubble center.
! Next:
!   - Direct OpenACC kernels for inner (i,j,k) loops.
!   - Outer bubble loop can remain sequential (small iteration count).
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
