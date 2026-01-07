# Kernel 351: s_uvw2rdr

## Source Location
- **File**: Src/uvw2rdr.f90
- **Subroutine**: s_uvw2rdr
- **Line**: ~166

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Apply analysis nudging forcing terms for velocity components

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: uvw2rdr.f90 :: s_uvw2rdr
! Summary : Apply analysis nudging forcing terms for velocity components
!           (u, v, w) to radar data with validity checks
! GPU diff: Medium
! Findings:
!   - Serial k-loop wrapping parallel i,j loops (private(k))
!   - Three separate conditional blocks for u, v, w components
!   - Inner conditional checks for valid radar data (lim34n threshold)
!   - Uses module variable lim34n from m_commath
! Next:
!   - Convert to OpenACC with collapse clause
!   - Inner conditionals may cause thread divergence on GPU
!   - Consider masking approach for better GPU efficiency
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
