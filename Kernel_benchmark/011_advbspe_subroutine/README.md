# Kernel 011: subroutine

## Source Location
- **File**: Src/advbspe.f90
- **Subroutine**: subroutine
- **Line**: ~125

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates base state pressure advection for horizontally

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: advbspe.f90 :: subroutine s_advbspe
! Summary : Calculates base state pressure advection for horizontally
!           explicit and vertically explicit method.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Uses module constant g from comphy.
!   - Pure arithmetic, all GPU compatible.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i).
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
