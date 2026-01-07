# Kernel 050: s_chkitr

## Source Location
- **File**: Src/chkitr.f90
- **Subroutine**: s_chkitr
- **Line**: ~144

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Find maximum absolute value of iteration variations (dvar)

## Runtime Profile (from test_real)
- **Calls**: 16
- **Average Loop Length**: 806.4K
- **Total Time**: 0.002s
- **Average Time per Call**: 0.133ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: chkitr.f90 :: s_chkitr
! Summary : Find maximum absolute value of iteration variations (dvar)
!           for convergence checking with max reduction
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function/subroutine calls inside parallel region
!   - Uses reduction(max:) on single variable intitc
!   - Uses intrinsic functions (abs, max)
!   - Simple 3D loop with element-wise max computation
!   - MPI_allreduce called after parallel region (not inside)
! Next:
!   - Direct OpenACC with collapse(3) and reduction(max:)
!   - GPU reduction primitives well-suited for this pattern
! Runtime:
!   - Calls: 16
!   - AvgLoops: 806.4K
!   - TotalTime: 0.002s (0.00%)
!   - AvgTime: 0.133ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
