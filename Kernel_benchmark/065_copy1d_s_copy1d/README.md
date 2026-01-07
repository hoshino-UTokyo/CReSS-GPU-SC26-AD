# Kernel 065: s_copy1d

## Source Location
- **File**: Src/copy1d.f90
- **Subroutine**: s_copy1d
- **Line**: ~101

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Simple 1D array copy from invar to outvar.

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 128
- **Total Time**: 0.000s
- **Average Time per Call**: 0.012ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: copy1d.f90 :: s_copy1d
! Summary : Simple 1D array copy from invar to outvar.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls
!   - Trivial memory copy operation
!   - Independent element-wise operations
! Next:
!   - Straightforward GPU port; consider using cudaMemcpy or similar
!   - For small arrays, overhead may exceed benefit of GPU execution
!   - May be better to keep data resident on GPU and avoid copy
! Runtime:
!   - Calls: 1
!   - AvgLoops: 128
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.012ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
