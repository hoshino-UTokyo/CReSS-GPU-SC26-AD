# Kernel 311: s_steptund

## Source Location
- **File**: Src/steptund.f90
- **Subroutine**: s_steptund
- **Line**: ~241
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Set boundary conditions and coefficient matrix for soil/sea

## Runtime Profile (from test_real)
- **Calls**: 18
- **Average Loop Length**: 806.4K
- **Total Time**: 0.093s
- **Average Time per Call**: 5.154ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: steptund.f90 :: s_steptund
! Summary : Set boundary conditions and coefficient matrix for soil/sea
!           temperature tridiagonal equation solver
! GPU diff: Medium
! Findings:
!   - Multiple conditional branches based on land type and sfcopt
!   - Nested k-loops with inner i,j loops parallelized
!   - Array updates depend on land use classification
!   - No function calls within parallel region
! Next:
!   - Consider collapsing k,j,i loops for better GPU occupancy
!   - Use data directives for tundp, tundf, rr, ss, tt arrays
! Runtime:
!   - Calls: 18
!   - AvgLoops: 806.4K
!   - TotalTime: 0.093s (0.00%)
!   - AvgTime: 5.154ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
