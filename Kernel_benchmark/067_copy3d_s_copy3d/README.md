# Kernel 067: s_copy3d

## Source Location
- **File**: Src/copy3d.f90
- **Subroutine**: s_copy3d
- **Line**: ~115

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Simple 3D array copy from invar to outvar.

## Runtime Profile (from test_real)
- **Calls**: 1446
- **Average Loop Length**: 103.9M
- **Total Time**: 3.996s
- **Average Time per Call**: 2.764ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: copy3d.f90 :: s_copy3d
! Summary : Simple 3D array copy from invar to outvar.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls
!   - Trivial memory copy operation
!   - Independent element-wise operations
!   - Outer k loop is sequential in current OpenMP structure
! Next:
!   - Collapse all three loops (k,j,i) for better GPU occupancy
!   - Consider using device-to-device memcpy for efficiency
!   - May be better to keep data resident on GPU and avoid copy calls
! Runtime:
!   - Calls: 1446
!   - AvgLoops: 103.9M
!   - TotalTime: 3.996s (0.13%)
!   - AvgTime: 2.764ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
