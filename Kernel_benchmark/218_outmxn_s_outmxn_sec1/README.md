# Kernel 218: s_outmxn

## Source Location
- **File**: Src/outmxn.f90
- **Subroutine**: s_outmxn
- **Line**: ~222
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Find max/min values of 3D variable using reductions with

## Runtime Profile (from test_real)
- **Calls**: 5415
- **Average Loop Length**: 102.5M
- **Total Time**: 5.493s
- **Average Time per Call**: 1.014ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: outmxn.f90 :: s_outmxn
! Summary : Find max/min values of 3D variable using reductions with
!           eps offset for numerical stability in comparisons.
! GPU diff: Medium
! Findings:
!   - Uses reduction(max/min) for maxvl, maxeps, minvl, mineps
!   - Triple nested loop over full 3D domain (i,j,k)
!   - Uses sign intrinsic for eps offset calculation
!   - No function calls; simple arithmetic operations
! Next:
!   - Data managed automatically via Unified Memory atomic or parallel reduction
!   - Consider using CUB or Thrust for reduction primitives
!   - May need two-pass approach for value then indices
! Runtime:
!   - Calls: 5415
!   - AvgLoops: 102.5M
!   - TotalTime: 5.493s (0.18%)
!   - AvgTime: 1.014ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
