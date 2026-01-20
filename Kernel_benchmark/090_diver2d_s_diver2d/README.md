# Kernel 090: s_diver2d

## Source Location
- **File**: Src/diver2d.f90
- **Subroutine**: s_diver2d
- **Line**: ~168

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate 2D horizontal negative divergence using velocity

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 102.5M
- **Total Time**: 162.333s
- **Average Time per Call**: 11.273ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diver2d.f90 :: s_diver2d
! Summary : Calculate 2D horizontal negative divergence using velocity
!           components u,v weighted by Jacobian and map scale factors.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - No global/module variable writes
!   - No synchronization constructs
!   - Multiple branches (mfcopt, mpopt) but all are simple data-parallel loops
!   - Two-phase computation: multiply then difference
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - tmp1, tmp2 are temporary arrays that can be fused or kept on GPU
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 102.5M
!   - TotalTime: 162.333s (5.45%)
!   - AvgTime: 11.273ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
