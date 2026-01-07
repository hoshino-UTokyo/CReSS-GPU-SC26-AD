# Kernel 077: s_defomten

## Source Location
- **File**: Src/defomten.f90
- **Subroutine**: s_defomten
- **Line**: ~342
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Calculate all components of deformation tensor (s11, s22, s33,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: defomten.f90 :: s_defomten (second parallel region)
! Summary : Calculate all components of deformation tensor (s11, s22, s33,
!           s12, s13, s23, s31, s32) using finite differences and Jacobians.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Private variable k for outer loop; jcbiv2, mfdvj2 for local scalars
!   - Writes to s11, s22, s33, s12, s13, s23, s31, s32, tmp1-tmp4 arrays
!   - Complex conditional branches based on trnopt, mfcopt, mpopt options
!   - Multiple sequential k-loops with data dependencies
!   - Stencil operations with varying grid indices
! Next:
!   - Split into multiple GPU kernels based on logical sections
!   - Create data region encompassing all arrays
!   - Consider kernel fusion where dependencies allow
!   - Profile to identify performance-critical sections
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
