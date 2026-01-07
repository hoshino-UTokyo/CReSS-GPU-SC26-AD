# Kernel 107: subroutine

## Source Location
- **File**: Src/exbcw.f90
- **Subroutine**: subroutine
- **Line**: ~211

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Forces lateral boundary values of w velocity to external

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: exbcw.f90 :: subroutine s_exbcw
! Summary : Forces lateral boundary values of w velocity to external
!           (GPV) boundary values with radiative/relaxation approach.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module variables (ebw, ebe, ebs, ebn, isub, jsub, etc.).
!   - No synchronization constructs.
!   - Uses intrinsic abs() - GPU compatible.
!   - Complex boundary-position-dependent conditionals (corners, edges).
!   - Only boundary cells are updated - sparse computation.
! Next:
!   - Consider separate kernels for each boundary region.
!   - Boundary-only work has low arithmetic intensity on GPU.
!   - Masking approach may be more efficient than conditionals.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
