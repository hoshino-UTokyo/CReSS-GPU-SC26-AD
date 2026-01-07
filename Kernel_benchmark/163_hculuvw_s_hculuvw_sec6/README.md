# Kernel 163: s_hculuvw

## Source Location
- **File**: Src/hculuvw.f90
- **Subroutine**: s_hculuvw
- **Line**: ~3414
- **Section**: 6 of 6 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Compute horizontal w-velocity advection using Cubic Lagrange scheme

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: hculuvw.f90 :: s_hculuvw
! Summary : Compute horizontal w-velocity advection using Cubic Lagrange scheme
!           with 4-point stencil and branch logic for wind direction.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Complex branching based on velocity sign (u8w, v8w directions)
!   - Reads from wp, writes to advd and wf arrays
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
