# Kernel 346: s_upwnbin

## Source Location
- **File**: Src/upwnbin.f90
- **Subroutine**: s_upwnbin
- **Line**: ~152

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate sedimentation for optional bin concentrations

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: upwnbin.f90 :: s_upwnbin
! Summary : Calculate sedimentation for optional bin concentrations
!           using upwind scheme with flux divergence
! GPU diff: Easy
! Findings:
!   - Serial k-loop wrapping parallel i,j loops (private(k))
!   - Three sequential loop nests: flux calc, update, boundary copy
!   - Contains max() intrinsic for non-negative concentration
!   - No complex branching or function calls
! Next:
!   - Collapse loops or use OpenACC kernels with loop directive
!   - Can potentially fuse kernels for better performance
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
