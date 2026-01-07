# Kernel 214: s_nuc2nd

## Source Location
- **File**: Src/nuc2nd.f90
- **Subroutine**: s_nuc2nd
- **Line**: ~140

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate secondary ice nucleation rate from snow and graupel

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 1.550s
- **Average Time per Call**: 0.034ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: nuc2nd.f90 :: s_nuc2nd
! Summary : Calculate secondary ice nucleation rate from snow and graupel
!           based on temperature and wet/dry graupel conditions.
! GPU diff: Easy
! Findings:
!   - Conditional branch for nk.eq.1 vs nk.gt.1 cases
!   - Temperature-based conditionals (270.16, 268.16, 265.16 K)
!   - Output arrays spsi, spgi written independently per grid point
!   - No function calls; simple arithmetic operations
!   - No inter-thread dependencies; fully parallel
! Next:
!   - Direct port to GPU kernel
!   - Branch divergence from temperature conditionals
!   - Consider using select case or predicated assignments
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.550s (0.05%)
!   - AvgTime: 0.034ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
