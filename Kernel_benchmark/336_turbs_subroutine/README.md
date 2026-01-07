# Kernel 336: subroutine

## Source Location
- **File**: Src/turbs.f90
- **Subroutine**: subroutine
- **Line**: ~193

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates scalar turbulent mixing (diffusion) with Jacobian

## Runtime Profile (from test_real)
- **Calls**: 3600
- **Average Loop Length**: 100.5M
- **Total Time**: 48.054s
- **Average Time per Call**: 13.348ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: turbs.f90 :: subroutine s_turbs
! Summary : Calculates scalar turbulent mixing (diffusion) with Jacobian
!           and map scale factor corrections for terrain-following coords.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Multiple code paths based on trnopt, mpopt, mfcopt options.
!   - Multi-stage: (1) compute tmp1,tmp2 flux components,
!     (2) optional tmp3 terrain correction, (3) divergence to sfrc.
!   - Uses 2D map scale factor arrays (mf, rmf, rmf8u, rmf8v).
! Next:
!   - Can use OpenACC kernels for each loop nest.
!   - Many conditional paths - consider unifying with flag-based selection.
!   - Temporary arrays already allocated.
! Runtime:
!   - Calls: 3600
!   - AvgLoops: 100.5M
!   - TotalTime: 48.054s (1.61%)
!   - AvgTime: 13.348ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
