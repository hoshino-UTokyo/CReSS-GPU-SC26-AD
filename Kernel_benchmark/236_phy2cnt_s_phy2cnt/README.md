# Kernel 236: s_phy2cnt

## Source Location
- **File**: Src/phy2cnt.f90
- **Subroutine**: s_phy2cnt
- **Line**: ~197

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate zeta components of contravariant velocity from

## Runtime Profile (from test_real)
- **Calls**: 15121
- **Average Loop Length**: 101.6M
- **Total Time**: 42.305s
- **Average Time per Call**: 2.798ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: phy2cnt.f90 :: s_phy2cnt
! Summary : Calculate zeta components of contravariant velocity from
!           physical velocity components with terrain-following coordinates.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to output array wc, work arrays j31u2, j32v2, mf25
!   - Simple 3D stencil computations with straightforward data access
!   - Multiple conditional paths based on trnopt, sthopt, mfcopt, mpopt
!   - No explicit barriers but implicit at !$omp end do
! Next:
!   - Use OpenACC parallel loop with collapse(3) for 3D loops
!   - Straightforward GPU port with data region for arrays
!   - Consider kernel fusion for consecutive loops
! Runtime:
!   - Calls: 15121
!   - AvgLoops: 101.6M
!   - TotalTime: 42.305s (1.42%)
!   - AvgTime: 2.798ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
