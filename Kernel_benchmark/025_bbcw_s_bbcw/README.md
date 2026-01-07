# Kernel 025: s_bbcw

## Source Location
- **File**: Src/bbcw.f90
- **Subroutine**: s_bbcw
- **Line**: ~144

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Set bottom boundary conditions for vertical velocity (wf)

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 807.3K
- **Total Time**: 1.749s
- **Average Time per Call**: 0.121ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bbcw.f90 :: s_bbcw
! Summary : Set bottom boundary conditions for vertical velocity (wf)
!           using terrain-following coordinate transformations
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to intent(inout) arrays wf, j31u2, j32v2
!   - Multiple conditional branches based on mfcopt and mpopt options
!   - Uses work arrays j31u2, j32v2 for intermediate calculations
! Next:
!   - GPU port requires handling conditional branches
!   - Consider separating compute kernels by mfcopt/mpopt case
!   - Work arrays j31u2, j32v2 should be device-resident
!   - Branch divergence from mfcopt/mpopt may impact performance
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 807.3K
!   - TotalTime: 1.749s (0.06%)
!   - AvgTime: 0.121ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
