# Kernel 013: subroutine

## Source Location
- **File**: Src/advbspt.f90
- **Subroutine**: subroutine
- **Line**: ~147

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates base state potential temperature advection by

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 101.2M
- **Total Time**: 3.017s
- **Average Time per Call**: 8.380ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: advbspt.f90 :: subroutine s_advbspt
! Summary : Calculates base state potential temperature advection by
!           vertical velocity for gravity wave mode calculations.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Two stages: (1) compute pta8w at w-points, (2) accumulate to ptadv.
!   - Conditional on gwmopt for accumulation vs assignment.
!   - All grid points are independent within each stage.
! Next:
!   - Direct OpenACC kernels for each loop nest.
!   - Consider fusing stages if pta8w is temporary.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 101.2M
!   - TotalTime: 3.017s (0.10%)
!   - AvgTime: 8.380ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
