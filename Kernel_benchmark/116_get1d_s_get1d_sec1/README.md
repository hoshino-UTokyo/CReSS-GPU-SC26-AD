# Kernel 116: s_get1d

## Source Location
- **File**: Src/get1d.f90
- **Subroutine**: s_get1d
- **Line**: ~293
- **Section**: 1 of 3 in this subroutine

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Compute horizontally averaged 1D profiles of velocity, pressure,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: get1d.f90 :: s_get1d
! Summary : Compute horizontally averaged 1D profiles of velocity, pressure,
!           temperature, and moisture by summing over grid and reducing.
! GPU diff: Hard
! Findings:
!   - Multiple reduction operations (+ for pref, u1d, v1d, pt1d, qv1d)
!   - Complex conditional logic based on fproc, idstr<=idend, gpvvar
!   - Grid point counting with ngcnt variable in some branches
!   - Many omp do regions with schedule(runtime)
!   - Accumulation into 1D arrays indexed by k_sub
! Next:
!   - Requires GPU reduction kernels for horizontal averaging
!   - May need multiple kernel launches for different conditional branches
!   - Consider restructuring to reduce conditional complexity
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
