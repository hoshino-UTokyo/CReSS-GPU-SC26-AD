# Kernel 086: s_diffequa

## Source Location
- **File**: Src/diffequa.f90
- **Subroutine**: s_diffequa
- **Line**: ~207

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Compute Laplacian terms for parabolic PDE iteration solving,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diffequa.f90 :: s_diffequa
! Summary : Compute Laplacian terms for parabolic PDE iteration solving,
!           updating phi with finite difference stencil in x, y, z.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - No global/module variable writes
!   - No synchronization constructs inside parallel region
!   - Simple 3D stencil computation with neighbor accesses
!   - Part of iterative solver (iterate loop outside) with MPI communication after parallel region
! Next:
!   - GPU offload possible but requires data management for iterative loop
!   - Consider keeping data on GPU across iterations to reduce transfer overhead
!   - MPI communication after parallel region needs attention for GPU-aware MPI
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
