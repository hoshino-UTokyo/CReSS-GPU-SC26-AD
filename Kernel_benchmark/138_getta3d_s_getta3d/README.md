# Kernel 138: s_getta3d

## Source Location
- **File**: Src/getta3d.f90
- **Subroutine**: s_getta3d
- **Line**: ~110

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute air temperature from base state potential temperature,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getta3d.f90 :: s_getta3d
! Summary : Compute air temperature from base state potential temperature,
!           perturbation, and Exner function at all grid points.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Simple element-wise arithmetic: t = (ptbr + ptp) * pi
!   - No global writes, only output array t is modified
!   - No synchronization constructs other than implicit barrier at end do
! Next:
!   - Direct translation to OpenACC or OpenACC with collapsed loops
!   - Consider loop collapse for k,j,i dimensions
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
