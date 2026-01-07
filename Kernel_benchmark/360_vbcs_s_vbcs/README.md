# Kernel 360: s_vbcs

## Source Location
- **File**: Src/vbcs.f90
- **Subroutine**: s_vbcs
- **Line**: ~116

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets vertical boundary conditions for scalar variable by copying

## Runtime Profile (from test_real)
- **Calls**: 3962
- **Average Loop Length**: 806.4K
- **Total Time**: 0.164s
- **Average Time per Call**: 0.042ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vbcs.f90 :: s_vbcs
! Summary : Sets vertical boundary conditions for scalar variable by copying
!           values from adjacent levels at bottom and top boundaries.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global/module variable writes, only local array writes
!   - No synchronization constructs (barrier, critical, atomic)
!   - Two separate omp do regions for bottom and top boundaries
! Next:
!   - Direct conversion to OpenACC parallel loop or OpenACC
!   - Both loops are independent and can run concurrently on GPU
! Runtime:
!   - Calls: 3962
!   - AvgLoops: 806.4K
!   - TotalTime: 0.164s (0.01%)
!   - AvgTime: 0.042ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
