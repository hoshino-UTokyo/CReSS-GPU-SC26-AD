# Kernel 054: s_chkrain

## Source Location
- **File**: Src/chkrain.f90
- **Subroutine**: s_chkrain
- **Line**: ~134

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Set precipitation flag (fall) based on water/ice precipitation

## Runtime Profile (from test_real)
- **Calls**: 361
- **Average Loop Length**: 806.4K
- **Total Time**: 0.025s
- **Average Time per Call**: 0.069ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: chkrain.f90 :: s_chkrain
! Summary : Set precipitation flag (fall) based on water/ice precipitation
!           thresholds for various microphysics options (bulk/bin methods)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function/subroutine calls inside parallel region
!   - No reductions or synchronization constructs
!   - Multiple conditional branches based on cphopt, haiopt, fmois flags
!   - Simple 2D loops setting output array fall to 1.0 or -1.0
!   - Reads from prwtr and price arrays
! Next:
!   - Direct OpenACC with collapse(2) for GPU
!   - Branching within kernel may cause thread divergence; consider separate kernels
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.025s (0.00%)
!   - AvgTime: 0.069ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
