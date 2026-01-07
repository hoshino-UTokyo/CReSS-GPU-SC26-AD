# Kernel 024: s_baserho

## Source Location
- **File**: Src/baserho.f90
- **Subroutine**: s_baserho
- **Line**: ~158

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Multiply base state density by Jacobian to compute rst array

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 102.9M
- **Total Time**: 0.004s
- **Average Time per Call**: 3.750ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: baserho.f90 :: s_baserho
! Summary : Multiply base state density by Jacobian to compute rst array
!           for use in atmospheric dynamics calculations
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to intent(out) array rst
!   - Simple element-wise computation: rst = abs(jcb) * rbr
!   - Outer k loop with inner parallel i,j loops
! Next:
!   - Straightforward GPU port with OpenACC parallel loops
!   - Collapse all three loops (k,j,i) for maximum parallelism
!   - Intrinsic abs function is GPU-compatible
! Runtime:
!   - Calls: 1
!   - AvgLoops: 102.9M
!   - TotalTime: 0.004s (0.00%)
!   - AvgTime: 3.750ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
