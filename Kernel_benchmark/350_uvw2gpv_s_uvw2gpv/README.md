# Kernel 350: s_uvw2gpv

## Source Location
- **File**: Src/uvw2gpv.f90
- **Subroutine**: s_uvw2gpv
- **Line**: ~183

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Apply analysis nudging forcing terms for velocity components

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: uvw2gpv.f90 :: s_uvw2gpv
! Summary : Apply analysis nudging forcing terms for velocity components
!           (u, v, w) to GPV data with time interpolation
! GPU diff: Easy
! Findings:
!   - Serial k-loop wrapping parallel i,j loops (private(k))
!   - Three separate conditional blocks for u, v, w components
!   - Conditionals based on character flag nggvar (checked outside loop)
!   - Simple accumulation to forcing arrays
! Next:
!   - Convert to OpenACC with collapse clause
!   - Conditionals are at outer level, no branch divergence in kernel
!   - Can separate into three independent kernels
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
