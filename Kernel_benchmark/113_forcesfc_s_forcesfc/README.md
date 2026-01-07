# Kernel 113: s_forcesfc

## Source Location
- **File**: Src/forcesfc.f90
- **Subroutine**: s_forcesfc
- **Line**: ~163

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute surface flux forcing terms for potential temperature,

## Runtime Profile (from test_real)
- **Calls**: 361
- **Average Loop Length**: 806.4K
- **Total Time**: 0.050s
- **Average Time per Call**: 0.139ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: forcesfc.f90 :: s_forcesfc
! Summary : Compute surface flux forcing terms for potential temperature,
!           water vapor, and velocity components at bottom boundary.
! GPU diff: Easy
! Findings:
!   - Multiple omp do regions for ptfrc, qvfrc, ufrc, vfrc calculations
!   - Conditional branches based on fmois (dry/moist) flag
!   - Uses sqrt intrinsic function for velocity calculations
!   - All operations are on 2D surface layer (k=1 or k=2)
!   - No reductions or synchronization
! Next:
!   - Port as 2D GPU kernels for surface layer
!   - Handle dry/moist branching with separate kernels or compile-time flag
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.050s (0.00%)
!   - AvgTime: 0.139ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
