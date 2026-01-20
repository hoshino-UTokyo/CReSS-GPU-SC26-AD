# Kernel 051: s_chkmoist

## Source Location
- **File**: Src/chkmoist.f90
- **Subroutine**: s_chkmoist
- **Line**: ~124

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Find maximum water vapor mixing ratio to determine if atmosphere

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 102.4M
- **Total Time**: 0.001s
- **Average Time per Call**: 1.117ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: chkmoist.f90 :: s_chkmoist
! Summary : Find maximum water vapor mixing ratio to determine if atmosphere
!           is moist or dry with max reduction
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function/subroutine calls inside parallel region
!   - Uses reduction(max:) on single variable qvmax
!   - Uses intrinsic function (max)
!   - Simple 3D loop with element-wise max computation
!   - MPI_allreduce called after parallel region (not inside)
! Next:
!   - Direct OpenACC with collapse(3) and reduction(max:)
!   - GPU reduction primitives well-suited for this pattern
! Runtime:
!   - Calls: 1
!   - AvgLoops: 102.4M
!   - TotalTime: 0.001s (0.00%)
!   - AvgTime: 1.117ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
