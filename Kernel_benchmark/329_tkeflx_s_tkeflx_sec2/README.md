# Kernel 329: s_tkeflx

## Source Location
- **File**: Src/tkeflx.f90
- **Subroutine**: s_tkeflx
- **Line**: ~282
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate x, y, z components of turbulent fluxes for TKE

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: tkeflx.f90 :: s_tkeflx (main flux calculation)
! Summary : Calculate x, y, z components of turbulent fluxes for TKE
!           with optional terrain and map scale factor corrections
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls within loops
!   - Writes to h1, h2, h3, jcbtke arrays
!   - Multiple conditional branches (mfcopt, mpopt, trnopt)
!   - h3 used as temporary array for rmf*rkh in some branches
!   - No synchronization constructs within parallel region
! Next:
!   - Multiple kernel approach based on options or unified kernel with conditionals
!   - Map scale factor arrays need to be on device
!   - Consider kernel specialization for common option combinations
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
