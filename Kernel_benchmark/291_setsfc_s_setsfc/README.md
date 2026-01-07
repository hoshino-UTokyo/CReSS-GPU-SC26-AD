# Kernel 291: s_setsfc

## Source Location
- **File**: Src/setsfc.f90
- **Subroutine**: s_setsfc
- **Line**: ~238

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates surface parameters including pressure, temperature, virtual potential

## Runtime Profile (from test_real)
- **Calls**: 361
- **Average Loop Length**: 102.4M
- **Total Time**: 4.230s
- **Average Time per Call**: 11.717ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: setsfc.f90 :: s_setsfc
! Summary : Calculates surface parameters including pressure, temperature, virtual potential
!           temperature, saturation mixing ratio, and velocity magnitude at lowest levels.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses intrinsic functions: exp, log, max, min, sqrt
!   - Uses module constants from m_commath and m_comphy (es0, t0, epsva, etc.)
!   - Multiple conditional branches based on fmois flag and land type
!   - Complex saturation vapor pressure calculations with exponentials
!   - All loops independent with private i,j,k and local scalar variables
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Port exp/log intrinsics directly (GPU-compatible)
!   - Consider separate kernels for dry vs moist branches
!   - Use data regions to minimize transfers of large 3D arrays
! Runtime:
!   - Calls: 361
!   - AvgLoops: 102.4M
!   - TotalTime: 4.230s (0.14%)
!   - AvgTime: 11.717ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
