# Kernel 005: subroutine

## Source Location
- **File**: Src/adjstnsg.f90
- **Subroutine**: subroutine
- **Line**: ~181

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Adjusts concentrations of ice precipitation (snow, graupel,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: adjstnsg.f90 :: subroutine s_adjstnsg
! Summary : Adjusts concentrations of ice precipitation (snow, graupel,
!           hail) to be within physical bounds using diagnostic formulas.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - Calls getiname() before parallel region (not inside).
!   - Uses module constants from comphy/commath.
!   - Uses intrinsic sqrt, min, max - all GPU compatible.
!   - Conditional on haiopt determines if hail is processed.
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
