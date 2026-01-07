# Kernel 144: s_getz

## Source Location
- **File**: Src/getz.f90
- **Subroutine**: s_getz
- **Line**: ~114

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate 1D zeta (terrain-following) vertical coordinates

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 128
- **Total Time**: 0.000s
- **Average Time per Call**: 0.016ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getz.f90 :: s_getz
! Summary : Calculate 1D zeta (terrain-following) vertical coordinates
!           from sea surface height and grid spacing.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic real function (GPU compatible)
!   - Simple 1D loop with arithmetic: z = zsfc + (k-2)*dz
!   - No global writes, only output array z is modified
! Next:
!   - 1D array with nk elements (typically small, <100)
!   - May not benefit from GPU offload due to small size
!   - If needed, use OpenACC with single team
! Runtime:
!   - Calls: 1
!   - AvgLoops: 128
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.016ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
