# Kernel 128: s_getmxn

## Source Location
- **File**: Src/getmxn.f90
- **Subroutine**: s_getmxn
- **Line**: ~257
- **Section**: 2 of 4 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Identifies grid indices (i,j,k) where the maximum value occurs

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getmxn.f90 :: s_getmxn
! Summary : Identifies grid indices (i,j,k) where the maximum value occurs
!           by comparing against the found maximum with tolerance
! GPU diff: Medium
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic abs(), max(), sign() functions - GPU compatible
!   - Contains reduction operations (max: maxi, maxj, maxk)
!   - Conditional update inside loop based on value comparison
!   - Small loop size (npe processors)
! Next:
!   - Small loop size may not benefit from GPU offloading
!   - Integer reductions supported in OpenACC
!   - Consider keeping on CPU due to small iteration count
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
