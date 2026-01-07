# Kernel 141: s_gettrn

## Source Location
- **File**: Src/gettrn.f90
- **Subroutine**: s_gettrn
- **Line**: ~263
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Generate bell-shaped mountain terrain using Gaussian-like

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: gettrn.f90 :: s_gettrn (trnopt=1 branch)
! Summary : Generate bell-shaped mountain terrain using Gaussian-like
!           formula with configurable height, width, and center position.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic max function (GPU compatible)
!   - Reads from 1D arrays xs(i), ys(j) - need to ensure GPU accessible
!   - Simple arithmetic with division and max
!   - No global writes, only output array ht is modified
! Next:
!   - Direct translation to OpenACC with teams distribute
!   - Ensure xs and ys arrays are mapped to device
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
