# Kernel 069: subroutine

## Source Location
- **File**: Src/coriuv.f90
- **Subroutine**: subroutine
- **Line**: ~126

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates Coriolis force terms for u and v velocity equations

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 100.6M
- **Total Time**: 5.402s
- **Average Time per Call**: 15.006ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: coriuv.f90 :: subroutine s_coriuv
! Summary : Calculates Coriolis force terms for u and v velocity equations
!           using f-plane or beta-plane approximation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Two-stage calculation: (1) compute tmp1, (2) add to forcing.
!   - Stencil averaging of velocity for Coriolis term.
!   - All grid points are independent within each loop nest.
! Next:
!   - Direct OpenACC kernels with collapse(2) on i,j loops.
!   - Consider fusing the two stages into single kernel per component.
!   - Temporary array tmp1 needed for staggered grid averaging.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.6M
!   - TotalTime: 5.402s (0.18%)
!   - AvgTime: 15.006ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
