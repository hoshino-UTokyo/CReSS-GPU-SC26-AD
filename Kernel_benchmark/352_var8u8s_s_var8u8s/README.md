# Kernel 352: s_var8u8s

## Source Location
- **File**: Src/var8u8s.f90
- **Subroutine**: s_var8u8s
- **Line**: ~103

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Average optional variable from u-points to scalar points

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: var8u8s.f90 :: s_var8u8s
! Summary : Average optional variable from u-points to scalar points
!           using simple 2-point stencil in x-direction
! GPU diff: Easy
! Findings:
!   - Serial k-loop wrapping parallel i,j loops (private(k))
!   - Simple 2-point averaging: (var8u(i,j,k)+var8u(i+1,j,k))*0.5
!   - No function calls or complex operations
!   - Output array is independent of input (no race condition)
! Next:
!   - Convert to OpenACC with collapse clause
!   - Straightforward GPU port with good memory access pattern
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
