# Kernel 106: s_exbcv

## Source Location
- **File**: Src/exbcv.f90
- **Subroutine**: s_exbcv
- **Line**: ~213

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Force lateral boundary values to external GPV boundary values

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 112.2K
- **Total Time**: 2.969s
- **Average Time per Call**: 0.206ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: exbcv.f90 :: s_exbcv
! Summary : Force lateral boundary values to external GPV boundary values
!           for y-component velocity using radiation boundary conditions
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Uses MPI domain decomposition variables (ebw, ebe, ebs, ebn, isub, jsub)
!   - Many conditional branches based on boundary location and options
!   - Processes south, north, west, east boundaries separately
!   - Updates v array at domain boundaries only
!   - Uses separate damping coefficients for tangential (exnews) and normal (exnorm)
!   - exbvar character flags control which boundaries are active
! Next:
!   - Boundary-only operations may not benefit much from GPU
!   - Consider keeping boundary conditions on CPU if main computation on GPU
!   - If porting, need separate small kernels for each boundary section
!   - MPI communication patterns need careful handling with GPU buffers
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 112.2K
!   - TotalTime: 2.969s (0.10%)
!   - AvgTime: 0.206ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
