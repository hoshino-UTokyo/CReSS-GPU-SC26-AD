# Kernel 014: subroutine

## Source Location
- **File**: Src/advp.f90
- **Subroutine**: subroutine
- **Line**: ~251

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates pressure advection using 2nd or 4th order finite

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.5M
- **Total Time**: 18.870s
- **Average Time per Call**: 52.416ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: advp.f90 :: subroutine s_advp
! Summary : Calculates pressure advection using 2nd or 4th order finite
!           difference methods with Jacobian weighting for terrain-following
!           coordinates.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region (pure arithmetic only).
!   - Reads module variable fourd3 (constant) from commath.
!   - No synchronization constructs.
!   - Multiple code paths based on advopt, mpopt, mfcopt, diaopt options.
!   - Stencil computations with temporary arrays (tmp1, tmp2, tmp3).
!   - Data dependency: 2nd order results feed into 4th order computation.
!   - All loops are embarrassingly parallel within each k-level.
! Next:
!   - Can use OpenACC parallel loop with collapse for (i,j) loops.
!   - Temporary arrays already allocated - good for GPU data management.
!   - Consider fusing kernels where possible to reduce memory traffic.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.5M
!   - TotalTime: 18.870s (0.63%)
!   - AvgTime: 52.416ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
