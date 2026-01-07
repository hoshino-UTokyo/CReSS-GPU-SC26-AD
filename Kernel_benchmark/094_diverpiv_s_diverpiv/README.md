# Kernel 094: s_diverpiv

## Source Location
- **File**: Src/diverpiv.f90
- **Subroutine**: s_diverpiv
- **Line**: ~122

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate vertical divergence for pressure equation (HEVI method),

## Runtime Profile (from test_real)
- **Calls**: 28800
- **Average Loop Length**: 100.4M
- **Total Time**: 105.301s
- **Average Time per Call**: 3.656ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diverpiv.f90 :: s_diverpiv
! Summary : Calculate vertical divergence for pressure equation (HEVI method),
!           computing rcsq * (w(k) - w(k+1)) * dziv at each grid point.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - No global/module variable writes
!   - No synchronization constructs
!   - Simple 3D loop with straightforward vertical differencing
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - Very simple kernel, good candidate for early GPU porting
! Runtime:
!   - Calls: 28800
!   - AvgLoops: 100.4M
!   - TotalTime: 105.301s (3.53%)
!   - AvgTime: 3.656ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
