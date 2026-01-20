# Kernel 002: subroutine

## Source Location
- **File**: Src/adjstnci.f90
- **Subroutine**: subroutine
- **Line**: ~123

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Adjusts cloud ice concentrations to be within physical bounds

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: adjstnci.f90 :: subroutine s_adjstnci
! Summary : Adjusts cloud ice concentrations to be within physical bounds
!           based on mixing ratio using min/max constraints.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Uses module constants (mimax, mi0) from comphy.
!   - Pure min/max operations, all GPU compatible.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i).
!   - Module constants passed as scalars to device.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
