# Kernel 045: subroutine

## Source Location
- **File**: Src/buoywse.f90
- **Subroutine**: subroutine
- **Line**: ~147

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates buoyancy forcing for small time step integration

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: buoywse.f90 :: subroutine s_buoywse
! Summary : Calculates buoyancy forcing for small time step integration
!           in horizontally/vertically explicit method.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constant g from comphy.
!   - No synchronization constructs.
!   - Conditional on gwmopt for buoyancy formulation selection.
!   - Two stages: (1) compute wb8s, (2) accumulate to wsml.
!   - All grid points are independent within each stage.
! Next:
!   - Direct OpenACC kernels for each loop nest.
!   - Consider fusing the two stages if wb8s is temporary.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
