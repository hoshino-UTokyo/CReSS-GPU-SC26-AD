# Kernel 079: s_depsitbw

## Source Location
- **File**: Src/depsitbw.f90
- **Subroutine**: s_depsitbw
- **Line**: ~288

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Perform bin-resolved deposition for water droplets, shifting

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: depsitbw.f90 :: s_depsitbw
! Summary : Perform bin-resolved deposition for water droplets, shifting
!           bin boundaries and updating mass/concentration distributions.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic exp, log, max, sqrt functions (GPU-compatible)
!   - Private variable n for bin category loop
!   - Writes to ptp, qv, mwbin, nwbin, bmws, mws, nws, ssw, lv, kp, dv, dm, etc.
!   - Complex bin microphysics with conditional logic per grid point
!   - Sequential dependencies: common vars -> bin shifts -> mass shifts -> remap
!   - Calls external subroutine remapbw after parallel region
! Next:
!   - Split into multiple GPU kernels: initialization, bin shifting, mass update
!   - Create persistent data region for bin arrays across kernels
!   - Handle remapbw call separately (may need separate GPU port)
!   - Profile to identify bottleneck loops
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
