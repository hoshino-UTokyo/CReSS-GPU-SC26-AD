# Kernel 349: subroutine

## Source Location
- **File**: Src/upwqp.f90
- **Subroutine**: subroutine
- **Line**: ~168

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates sedimentation flux and precipitation for optional

## Runtime Profile (from test_real)
- **Calls**: 1800
- **Average Loop Length**: 102.4M
- **Total Time**: 17.374s
- **Average Time per Call**: 9.652ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: upwqp.f90 :: subroutine s_upwqp
! Summary : Calculates sedimentation flux and precipitation for optional
!           precipitation mixing ratio using upwind scheme.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region (pure arithmetic only).
!   - No writes to global/module variables.
!   - No synchronization constructs.
!   - Uses intrinsic max() which is GPU-compatible.
!   - Vertical dependency: qpflx computed first, then used for qpf update.
!   - precip accumulation has no race (each (i,j) independent).
! Next:
!   - Split into two kernels: (1) compute qpflx, (2) update qpf and precip.
!   - Or use OpenACC with proper data clauses.
! Runtime:
!   - Calls: 1800
!   - AvgLoops: 102.4M
!   - TotalTime: 17.374s (0.58%)
!   - AvgTime: 9.652ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
