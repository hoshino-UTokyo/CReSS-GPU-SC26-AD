# Kernel 385: subroutine

## Source Location
- **File**: Src/vspuvw.f90
- **Subroutine**: subroutine
- **Line**: ~214

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates vertical sponge damping for u, v, w velocity

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 101.6M
- **Total Time**: 8.704s
- **Average Time per Call**: 24.177ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vspuvw.f90 :: subroutine s_vspuvw
! Summary : Calculates vertical sponge damping for u, v, w velocity
!           components near model top to absorb outgoing waves.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Multiple code paths based on vspvar and vspopt flags.
!   - Loops start from ksp0 (sponge layer index) to top.
!   - All grid points are independent within each loop.
!   - Damping to GPV data or base state based on vspopt.
! Next:
!   - Direct OpenACC kernels for each component.
!   - Evaluate conditions outside kernel to select code path.
!   - Consider fusing loops for same component if beneficial.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 101.6M
!   - TotalTime: 8.704s (0.29%)
!   - AvgTime: 24.177ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
