# Kernel 211: s_nlsmuvw

## Source Location
- **File**: Src/nlsmuvw.f90
- **Subroutine**: s_nlsmuvw
- **Line**: ~185

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Non-linear smoothing for u, v, w velocity components using

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: nlsmuvw.f90 :: s_nlsmuvw
! Summary : Non-linear smoothing for u, v, w velocity components using
!           finite differences with a*abs(a) nonlinear diffusion.
! GPU diff: Medium
! Findings:
!   - Multiple sequential do-k loops for u, v, w components
!   - tmp4 reused for each velocity component (u, v, w)
!   - tmp1, tmp2, tmp3 store intermediate differences per component
!   - Separate forcing updates for ufrc, vfrc, wfrc
!   - No function calls; uses intrinsic abs only
!   - Complex stencil patterns with different index ranges per loop
! Next:
!   - Consider separate kernels for u, v, w smoothing
!   - May need 3+ kernel launches per velocity component
!   - Ensure tmp4 synchronization between velocity components
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
