# Kernel 053: s_chkmxn

## Source Location
- **File**: Src/chkmxn.f90
- **Subroutine**: s_chkmxn
- **Line**: ~284
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Find indices (i,j,k) of maximum and minimum values in 3D array

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: chkmxn.f90 :: s_chkmxn
! Summary : Find indices (i,j,k) of maximum and minimum values in 3D array
!           using max/min reductions on index variables
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function/subroutine calls inside parallel region
!   - Uses multiple reductions on indices: reduction(max:maxi,maxj,maxk) reduction(min:mini,minj,mink)
!   - Uses intrinsic functions (abs, sign, max, min)
!   - Depends on maxeps/mineps computed in previous parallel region
!   - Conditional processing based on fproc flag and undefined value range
! Next:
!   - Index-finding reductions can be tricky on GPU; consider argmax/argmin patterns
!   - May need custom reduction or atomic compare-and-swap for indices
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
