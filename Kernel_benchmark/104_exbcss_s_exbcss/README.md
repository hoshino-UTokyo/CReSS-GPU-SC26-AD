# Kernel 104: s_exbcss

## Source Location
- **File**: Src/exbcss.f90
- **Subroutine**: s_exbcss
- **Line**: ~215

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Force lateral boundary values to external GPV boundary values

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 125
- **Total Time**: 3.244s
- **Average Time per Call**: 0.225ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: exbcss.f90 :: s_exbcss
! Summary : Force lateral boundary values to external GPV boundary values
!           for scalar variables using radiation boundary conditions
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Uses MPI domain decomposition variables (ebw, ebe, ebs, ebn, isub, jsub)
!   - Many conditional branches based on boundary location and options
!   - Processes corners, west, east, south, north boundaries separately
!   - Updates s array at domain boundaries only
!   - Small time step integration (dts) for acoustic mode
!   - exbvar character flags control which boundaries are active
! Next:
!   - Boundary-only operations may not benefit much from GPU
!   - Consider keeping boundary conditions on CPU if main computation on GPU
!   - If porting, need separate small kernels for each boundary section
!   - MPI communication patterns need careful handling with GPU buffers
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 125
!   - TotalTime: 3.244s (0.11%)
!   - AvgTime: 0.225ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
