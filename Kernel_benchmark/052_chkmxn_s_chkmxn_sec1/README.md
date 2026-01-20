# Kernel 052: s_chkmxn

## Source Location
- **File**: Src/chkmxn.f90
- **Subroutine**: s_chkmxn
- **Line**: ~189
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Find maximum and minimum values of 3D data array with optional

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: chkmxn.f90 :: s_chkmxn
! Summary : Find maximum and minimum values of 3D data array with optional
!           undefined value filtering, using max/min reductions
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function/subroutine calls inside parallel region
!   - Uses multiple reductions: reduction(max:maxvl,maxeps) reduction(min:minvl,mineps)
!   - Uses intrinsic functions (sign, max, min)
!   - Conditional processing based on fproc flag and undefined value range
!   - Simple 3D loop with element-wise min/max computation
! Next:
!   - Direct OpenACC with collapse(3) and multiple reductions
!   - GPU reduction primitives well-suited for this pattern
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
