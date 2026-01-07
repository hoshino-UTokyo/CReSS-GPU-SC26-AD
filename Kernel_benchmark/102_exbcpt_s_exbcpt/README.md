# Kernel 102: s_exbcpt

## Source Location
- **File**: Src/exbcpt.f90
- **Subroutine**: s_exbcpt
- **Line**: ~241

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Force lateral boundary values to external GPV boundary values

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 1
- **Total Time**: 0.124s
- **Average Time per Call**: 0.345ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: exbcpt.f90 :: s_exbcpt
! Summary : Force lateral boundary values to external GPV boundary values
!           for potential temperature using radiation boundary conditions
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Uses MPI domain decomposition variables (ebw, ebe, ebs, ebn, isub, jsub)
!   - Many conditional branches based on boundary location and options
!   - Processes corners, west, east, south, north boundaries separately
!   - Updates ptpf array at domain boundaries only
!   - advopt controls time stepping scheme (leapfrog vs other)
!   - exbvar character flags control which boundaries are active
! Next:
!   - Boundary-only operations may not benefit much from GPU
!   - Consider keeping boundary conditions on CPU if main computation on GPU
!   - If porting, need separate small kernels for each boundary section
!   - MPI communication patterns need careful handling with GPU buffers
! Runtime:
!   - Calls: 360
!   - AvgLoops: 1
!   - TotalTime: 0.124s (0.00%)
!   - AvgTime: 0.345ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
