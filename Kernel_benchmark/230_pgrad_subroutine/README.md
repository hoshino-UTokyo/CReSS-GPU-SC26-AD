# Kernel 230: subroutine

## Source Location
- **File**: Src/pgrad.f90
- **Subroutine**: subroutine
- **Line**: ~278

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates pressure gradient force for u, v, w equations

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 102.4M
- **Total Time**: 292.883s
- **Average Time per Call**: 20.339ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: pgrad.f90 :: subroutine s_pgrad
! Summary : Calculates pressure gradient force for u, v, w equations
!           including divergence damping and terrain-following corrections.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - Calls diver3d() before parallel region - need to check that subroutine.
!   - Reads module constant divndc from comphy (indirectly via divch/divcv).
!   - No synchronization constructs within parallel region.
!   - Multiple code paths based on divopt, trnopt, mfcopt, mpopt.
!   - Multi-stage stencil with temporary arrays (tmp1, tmp2, tmp3).
!   - Terrain correction (trnopt>=1) adds extra stencil stage.
! Next:
!   - Ensure diver3d is GPU-ready before porting this routine.
!   - Use multiple kernels matching the loop structure.
!   - Map scale factor conditionals can be evaluated outside kernel.
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 102.4M
!   - TotalTime: 292.883s (9.83%)
!   - AvgTime: 20.339ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
