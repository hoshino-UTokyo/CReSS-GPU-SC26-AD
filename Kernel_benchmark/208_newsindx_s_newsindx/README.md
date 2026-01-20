# Kernel 208: s_newsindx

## Source Location
- **File**: Src/newsindx.f90
- **Subroutine**: s_newsindx
- **Line**: ~181

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Find min/max data indices for base state variables using

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: newsindx.f90 :: s_newsindx
! Summary : Find min/max data indices for base state variables using
!           floor/int operations on real indices ri, rj arrays.
! GPU diff: Easy
! Findings:
!   - Uses reduction(min/max) for cidstr/cidend/cjdstr/cjdend
!   - No function calls inside parallel region (floor/int are intrinsics)
!   - No global writes, only local variable updates
!   - Conditional branch (mpopt.lt.10) outside omp do regions
! Next:
!   - Data managed automatically via Unified Memory atomic or warp-level reductions
!   - Consider fusing the two branches into one kernel with masking
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
