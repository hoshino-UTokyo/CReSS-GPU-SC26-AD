# Kernel 087: s_disptke

## Source Location
- **File**: Src/disptke.f90
- **Subroutine**: s_disptke
- **Line**: ~181

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate TKE dissipation term using turbulent length scale,

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 100.4M
- **Total Time**: 2.484s
- **Average Time per Call**: 6.899ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: disptke.f90 :: s_disptke
! Summary : Calculate TKE dissipation term using turbulent length scale,
!           with different formulations for isotropic/anisotropic cases.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region (only intrinsics: exp, log, sqrt)
!   - No global/module variable writes
!   - No synchronization constructs
!   - Multiple branches (isoopt, mfcopt, mpopt) but all are data-parallel loops
!   - Accumulation into tkefrc (inout), but each (i,j,k) is independent
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - Branch conditions can be hoisted outside target region for cleaner GPU code
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.4M
!   - TotalTime: 2.484s (0.08%)
!   - AvgTime: 6.899ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
