# Kernel 253: subroutine

## Source Location
- **File**: Src/rbcss.f90
- **Subroutine**: subroutine
- **Line**: ~242

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Sets radiative lateral boundary conditions for optional scalar

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rbcss.f90 :: subroutine s_rbcss
! Summary : Sets radiative lateral boundary conditions for optional scalar
!           variables with phase speed advection and damping.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region (pure arithmetic only).
!   - Reads module variables (ebw, ebe, ebs, ebn, isub, jsub, etc.).
!   - No synchronization constructs.
!   - Complex conditional structure based on boundary position (corners,
!     edges) - may cause thread divergence on GPU.
!   - Each boundary region is independent but sparse (only boundary cells).
! Next:
!   - Consider separate kernels for each boundary region to reduce divergence.
!   - Alternatively, use masking approach for unified kernel.
!   - Boundary-only computation means low arithmetic intensity.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
