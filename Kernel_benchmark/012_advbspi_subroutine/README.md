# Kernel 012: subroutine

## Source Location
- **File**: Src/advbspi.f90
- **Subroutine**: subroutine
- **Line**: ~154

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates base state pressure advection for horizontally

## Runtime Profile (from test_real)
- **Calls**: 28800
- **Average Loop Length**: 100.4M
- **Total Time**: 172.825s
- **Average Time per Call**: 6.001ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: advbspi.f90 :: subroutine s_advbspi
! Summary : Calculates base state pressure advection for horizontally
!           explicit and vertically implicit method (forward/backward).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - Calls getrname() before parallel region (not inside).
!   - Uses module constant g from comphy.
!   - Conditional on fproc determines backward vs forward calculation.
!   - Pure arithmetic, all GPU compatible.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i).
!   - May split into two kernels for back/fore branches.
! Runtime:
!   - Calls: 28800
!   - AvgLoops: 100.4M
!   - TotalTime: 172.825s (5.80%)
!   - AvgTime: 6.001ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
