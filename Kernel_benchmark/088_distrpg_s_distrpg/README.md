# Kernel 088: s_distrpg

## Source Location
- **File**: Src/distrpg.f90
- **Subroutine**: s_distrpg
- **Line**: ~162

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Distribute collision rates between rain and snow to graupel,

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 1.393s
- **Average Time per Call**: 0.030ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: distrpg.f90 :: s_distrpg
! Summary : Distribute collision rates between rain and snow to graupel,
!           based on diameter ratios and temperature thresholds.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region (only intrinsic: abs)
!   - No global/module variable writes
!   - No synchronization constructs
!   - Multiple branches (nk, cphopt) but all loops are data-parallel
!   - Conditional updates per grid point (temperature check, threshold)
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - Conditionals inside loop are fine for GPU (divergent but manageable)
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.393s (0.05%)
!   - AvgTime: 0.030ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
