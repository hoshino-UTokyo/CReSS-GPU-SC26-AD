# Kernel 381: s_vspp

## Source Location
- **File**: Src/vspp.f90
- **Subroutine**: s_vspp
- **Line**: ~159

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Applies vertical sponge damping to pressure forcing term,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vspp.f90 :: s_vspp
! Summary : Applies vertical sponge damping to pressure forcing term,
!           either relaxing to GPV data or base state value.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Conditional branch (vspopt, gpvvar) selects damping target
!   - k loop starts from ksp0-1 (variable start index)
!   - Simple arithmetic update to pfrc array
!   - No synchronization constructs other than implicit barriers
! Next:
!   - Straightforward GPU port with collapse on j,i loops
!   - Handle variable k-range start with appropriate kernel bounds
!   - Map pfrc, ppp, ppgpv, pptd, rbct, jcb arrays to device
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
