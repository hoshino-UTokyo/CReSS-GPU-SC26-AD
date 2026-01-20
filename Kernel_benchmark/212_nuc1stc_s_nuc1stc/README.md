# Kernel 212: s_nuc1stc

## Source Location
- **File**: Src/nuc1stc.f90
- **Subroutine**: s_nuc1stc
- **Line**: ~194

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate ice nucleation rates (condensation, contact,

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 1.366s
- **Average Time per Call**: 0.030ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: nuc1stc.f90 :: s_nuc1stc
! Summary : Calculate ice nucleation rates (condensation, contact,
!           homogeneous) based on temperature and cloud water.
! GPU diff: Medium
! Findings:
!   - Conditional branch for nk.eq.1 vs nk.gt.1 cases
!   - Complex conditionals (temperature thresholds) per grid point
!   - Uses exp, log, min intrinsics - GPU compatible
!   - Many private variables (tc, piv, knd, dar, f1, f2, ft, etc.)
!   - Output array nuci written independently per grid point
!   - No inter-thread dependencies within each omp do
! Next:
!   - Branch divergence may reduce GPU efficiency
!   - Consider precomputing masks for temperature conditions
!   - Can port as single kernel with good occupancy
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.366s (0.05%)
!   - AvgTime: 0.030ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
