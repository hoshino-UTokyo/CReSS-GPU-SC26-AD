# Kernel 321: s_swp2nxt

## Source Location
- **File**: Src/swp2nxt.f90
- **Subroutine**: s_swp2nxt
- **Line**: ~355

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Swap prognostic variables (velocity, pressure, temperature,

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.5M
- **Total Time**: 30.923s
- **Average Time per Call**: 85.898ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: swp2nxt.f90 :: s_swp2nxt
! Summary : Swap prognostic variables (velocity, pressure, temperature,
!           hydrometeors, aerosols, tracers, TKE) between time levels
! GPU diff: Medium
! Findings:
!   - Many conditional branches based on advopt, cphopt, haiopt, etc.
!   - Simple array copy operations within loops
!   - Large number of arrays to swap (u,v,w,pp,ptp,qv,qwtr,qice,etc.)
!   - No function calls within parallel region
!   - Different swap patterns for centered vs Lagrange advection
! Next:
!   - Consider batching array swaps for GPU memory efficiency
!   - Use async data transfers if arrays already on GPU
!   - Collapse loops where possible for better occupancy
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.5M
!   - TotalTime: 30.923s (1.04%)
!   - AvgTime: 85.898ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
