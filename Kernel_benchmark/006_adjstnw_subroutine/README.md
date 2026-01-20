# Kernel 006: subroutine

## Source Location
- **File**: Src/adjstnw.f90
- **Subroutine**: subroutine
- **Line**: ~147

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Adjusts concentrations of cloud water and rain water to be

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: adjstnw.f90 :: subroutine s_adjstnw
! Summary : Adjusts concentrations of cloud water and rain water to be
!           within physical bounds using diagnostic relationships.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Uses module constants from comphy/commath.
!   - Uses intrinsic sqrt, min, max - all GPU compatible.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i).
!   - Module constants can be passed as scalars to device.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
