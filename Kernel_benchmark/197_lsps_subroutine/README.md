# Kernel 197: subroutine

## Source Location
- **File**: Src/lsps.f90
- **Subroutine**: subroutine
- **Line**: ~186

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Applies lateral sponge damping for scalar variables near

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: lsps.f90 :: subroutine s_lsps
! Summary : Applies lateral sponge damping for scalar variables near
!           domain boundaries, relaxing toward GPV or base state values.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Uses intrinsic mod() for option check (outside parallel).
!   - Two stages: (1) compute tmp1 (deviation), (2) apply damping.
!   - Optional smoothing stencil when lspopt >= 10.
!   - All grid points are independent within each stage.
! Next:
!   - Direct OpenACC kernels for each loop nest.
!   - rbcxy is 2D, can be efficiently accessed on GPU.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
