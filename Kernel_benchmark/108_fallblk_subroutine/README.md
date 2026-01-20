# Kernel 108: subroutine

## Source Location
- **File**: Src/fallblk.f90
- **Subroutine**: subroutine
- **Line**: ~315

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates minimum time step for precipitation fallout based on

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.4M
- **Total Time**: 2.581s
- **Average Time per Call**: 7.169ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: fallblk.f90 :: subroutine s_fallblk
! Summary : Calculates minimum time step for precipitation fallout based on
!           CFL condition with terminal velocities of hydrometeors.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constant eps from commath via temparam.
!   - Uses OpenMP reduction(min:) for dtpc, dtpr, dtpi, dtps, dtpg, dtph.
!   - Uses intrinsic min() - GPU compatible.
!   - Multiple code paths based on flqcqi_opt, haiopt.
!   - Note: After parallel region, calls upwqp, upwnp, upwqcg in a loop.
! Next:
!   - Reduction requires GPU-compatible reduction pattern.
!   - Use atomicMin or warp-level reduction for CFL check.
!   - Consider computing reductions in separate kernel before main loop.
!   - The upwqp/upwnp/upwqcg calls need separate GPU porting.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 2.581s (0.09%)
!   - AvgTime: 7.169ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
