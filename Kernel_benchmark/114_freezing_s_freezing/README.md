# Kernel 114: s_freezing

## Source Location
- **File**: Src/freezing.f90
- **Subroutine**: s_freezing
- **Line**: ~160

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate freezing rate from rain water to graupel based on

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 1.296s
- **Average Time per Call**: 0.028ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: freezing.f90 :: s_freezing
! Summary : Calculate freezing rate from rain water to graupel based on
!           temperature and rain water content with different cphopt modes.
! GPU diff: Medium
! Findings:
!   - Complex conditional logic (nk==1 vs nk>1, cphopt==2 vs >=3)
!   - Uses exp intrinsic function
!   - Multiple threshold checks (thresq, t0cel, tclow)
!   - Writes to frrg and frrgn arrays
!   - No reductions or synchronization
! Next:
!   - Port with GPU kernels handling conditionals via masks or separate kernels
!   - Consider separating cphopt==2 and cphopt>=3 cases for clarity
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.296s (0.04%)
!   - AvgTime: 0.028ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
