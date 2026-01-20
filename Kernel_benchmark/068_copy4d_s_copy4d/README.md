# Kernel 068: s_copy4d

## Source Location
- **File**: Src/copy4d.f90
- **Subroutine**: s_copy4d
- **Line**: ~125

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Simple 4D array copy from invar to outvar.

## Runtime Profile (from test_real)
- **Calls**: 3
- **Average Loop Length**: 277.1M
- **Total Time**: 0.026s
- **Average Time per Call**: 8.662ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: copy4d.f90 :: s_copy4d
! Summary : Simple 4D array copy from invar to outvar.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls
!   - Trivial memory copy operation
!   - Independent element-wise operations
!   - Outer n,k loops are sequential in current OpenMP structure
! Next:
!   - Collapse all four loops (n,k,j,i) for better GPU occupancy
!   - Consider using device-to-device memcpy for efficiency
!   - May be better to keep data resident on GPU and avoid copy calls
! Runtime:
!   - Calls: 3
!   - AvgLoops: 277.1M
!   - TotalTime: 0.026s (0.00%)
!   - AvgTime: 8.662ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
