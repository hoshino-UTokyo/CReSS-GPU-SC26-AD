# Kernel 266: s_rstuvwc

## Source Location
- **File**: Src/rstuvwc.f90
- **Subroutine**: s_rstuvwc
- **Line**: ~176

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Multiply base state density x Jacobian by velocity components u, v, and wc

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.5M
- **Total Time**: 4.050s
- **Average Time per Call**: 11.251ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rstuvwc.f90 :: s_rstuvwc
! Summary : Multiply base state density x Jacobian by velocity components u, v, and wc
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Multiple conditional branches based on mfcopt and mpopt
!   - Outer k loop is serial with inner !$omp do on i,j
!   - Three separate loop nests for rstxu, rstxv, rstxwc
!   - Simple element-wise multiplication and assignment
!   - No synchronization constructs
! Next:
!   - Collapse k,j,i loops for better GPU parallelism
!   - Consider separate kernels for u, v, wc computations
!   - Use OpenACC teams distribute parallel for collapse(3)
!   - Fuse the three loop nests if possible for better memory access
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.5M
!   - TotalTime: 4.050s (0.14%)
!   - AvgTime: 11.251ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
