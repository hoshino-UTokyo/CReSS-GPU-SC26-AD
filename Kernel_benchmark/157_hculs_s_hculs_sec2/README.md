# Kernel 157: s_hculs

## Source Location
- **File**: Src/hculs.f90
- **Subroutine**: s_hculs
- **Line**: ~294
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Compute horizontal scalar advection using Cubic Lagrange scheme

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: hculs.f90 :: s_hculs
! Summary : Compute horizontal scalar advection using Cubic Lagrange scheme
!           with 4-point stencil and branch logic for wind direction.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Complex branching based on velocity sign (u8s, v8s directions)
!   - Reads from sp, writes to advx and sf arrays
!   - Multiple code paths for mfcopt (map scale factor) options
!   - No sync constructs
! Next:
!   - Convert to OpenACC with collapse(2) on j,i loops
!   - Consider predicated execution or warp divergence mitigation
!   - Branch logic may cause GPU thread divergence
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
