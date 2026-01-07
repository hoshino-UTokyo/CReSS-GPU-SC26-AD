# Kernel 353: s_var8uvw

## Source Location
- **File**: Src/var8uvw.f90
- **Subroutine**: s_var8uvw
- **Line**: ~170

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Averages optional variable to u, v, and w staggered grid points

## Runtime Profile (from test_real)
- **Calls**: 2
- **Average Loop Length**: 102.8M
- **Total Time**: 0.016s
- **Average Time per Call**: 8.103ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: var8uvw.f90 :: s_var8uvw
! Summary : Averages optional variable to u, v, and w staggered grid points
!           using simple 2-point averaging in each direction.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global/module variable writes, only local array writes
!   - No synchronization constructs (barrier, critical, atomic)
!   - Simple loop structure with private loop indices
! Next:
!   - Direct conversion to OpenACC parallel loop or OpenACC
!   - Consider collapsing nested loops for better GPU utilization
! Runtime:
!   - Calls: 2
!   - AvgLoops: 102.8M
!   - TotalTime: 0.016s (0.00%)
!   - AvgTime: 8.103ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
