# Kernel 320: s_swadjst

## Source Location
- **File**: Src/swadjst.f90
- **Subroutine**: s_swadjst
- **Line**: ~237

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Perform saturation adjustment for water phase, updating

## Runtime Profile (from test_real)
- **Calls**: 720
- **Average Loop Length**: 102.4M
- **Total Time**: 13.608s
- **Average Time per Call**: 18.900ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: swadjst.f90 :: s_swadjst
! Summary : Perform saturation adjustment for water phase, updating
!           temperature, vapor, and cloud water with condensation/evaporation
! GPU diff: Hard
! Findings:
!   - Complex conditionals based on cphopt and saturation state
!   - Multiple exp/log intrinsic function calls per grid point
!   - Iterative two-pass adjustment for accuracy
!   - Deep nesting with many local temporary variables
!   - Conditional updates to ncc (cloud concentration) based on vertical velocity
! Next:
!   - Significant thread divergence expected from conditionals
!   - Consider separating cphopt<=3 and cphopt==4 paths
!   - Profile exp/log operations for GPU performance
! Runtime:
!   - Calls: 720
!   - AvgLoops: 102.4M
!   - TotalTime: 13.608s (0.46%)
!   - AvgTime: 18.900ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
