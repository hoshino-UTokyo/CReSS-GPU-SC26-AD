# Kernel 007: subroutine

## Source Location
- **File**: Src/adjstq.f90
- **Subroutine**: subroutine
- **Line**: ~137

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Forces hydrometeor mixing ratios (qv, qwtr, qice) to be

## Runtime Profile (from test_real)
- **Calls**: 3
- **Average Loop Length**: 102.4M
- **Total Time**: 0.043s
- **Average Time per Call**: 14.438ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: adjstq.f90 :: subroutine s_adjstq
! Summary : Forces hydrometeor mixing ratios (qv, qwtr, qice) to be
!           non-negative using max(0) for various cloud physics options.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - Calls getiname() before parallel region (not inside).
!   - Pure max() operations, all GPU compatible.
!   - Conditionals on cphopt/haiopt control which arrays are processed.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i).
!   - May split into separate kernels for different cphopt branches.
! Runtime:
!   - Calls: 3
!   - AvgLoops: 102.4M
!   - TotalTime: 0.043s (0.00%)
!   - AvgTime: 14.438ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
