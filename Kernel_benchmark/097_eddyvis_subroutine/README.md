# Kernel 097: subroutine

## Source Location
- **File**: Src/eddyvis.f90
- **Subroutine**: subroutine
- **Line**: ~252

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates eddy viscosity and turbulent length scale using

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 127
- **Total Time**: 5.280s
- **Average Time per Call**: 14.667ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: eddyvis.f90 :: subroutine s_eddyvis
! Summary : Calculates eddy viscosity and turbulent length scale using
!           Smagorinsky or Deardorff (TKE-based) formulations with
!           stability corrections.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constants (oned3, csnum, ckm, ckmin, prnum, kappa, eps)
!     from commath/comphy.
!   - No synchronization constructs.
!   - Uses intrinsic abs(), exp(), log(), max(), min(), sqrt() - all GPU ok.
!   - Complex nested conditionals (tubopt, isoopt, sfcopt, mfcopt, mpopt).
!   - Thread divergence likely due to many branching paths.
!   - All grid points are independent within selected code path.
! Next:
!   - Consider restructuring conditionals for GPU - evaluate outside kernel.
!   - Use template/variant approach for different physics configurations.
!   - Intrinsic functions are well-supported on GPU.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 127
!   - TotalTime: 5.280s (0.18%)
!   - AvgTime: 14.667ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
