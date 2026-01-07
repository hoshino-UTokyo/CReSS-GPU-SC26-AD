# Kernel 363: s_vbcw

## Source Location
- **File**: Src/vbcw.f90
- **Subroutine**: s_vbcw
- **Line**: ~184

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Sets vertical boundary conditions for z-velocity component using

## Runtime Profile (from test_real)
- **Calls**: 14401
- **Average Loop Length**: 806.4K
- **Total Time**: 5.081s
- **Average Time per Call**: 0.353ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vbcw.f90 :: s_vbcw
! Summary : Sets vertical boundary conditions for z-velocity component using
!           terrain-following coordinate transformations at bottom/top.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global/module variable writes, only local array writes
!   - No synchronization constructs (barrier, critical, atomic)
!   - Complex conditional branching on bbc, tbc, mfcopt, mpopt values
!   - Many separate omp do regions (20+) within single parallel region
!   - Uses temporary 2D arrays (mf25, j31u2, j32v2) for intermediate results
!   - Data dependency: j31u2/j32v2 computed then used in subsequent loops
! Next:
!   - Consider restructuring to reduce number of kernel launches on GPU
!   - Ensure proper data movement for intermediate 2D arrays
!   - May benefit from fusing some loops where data dependencies allow
!   - Conditionals can be handled with masked operations or separate kernels
! Runtime:
!   - Calls: 14401
!   - AvgLoops: 806.4K
!   - TotalTime: 5.081s (0.17%)
!   - AvgTime: 0.353ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
