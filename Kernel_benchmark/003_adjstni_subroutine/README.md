# Kernel 003: subroutine

## Source Location
- **File**: Src/adjstni.f90
- **Subroutine**: subroutine
- **Line**: ~192

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Adjusts concentrations of ice hydrometeors (cloud ice, snow,

## Runtime Profile (from test_real)
- **Calls**: 1080
- **Average Loop Length**: 102.4M
- **Total Time**: 13.236s
- **Average Time per Call**: 12.255ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: adjstni.f90 :: subroutine s_adjstni
! Summary : Adjusts concentrations of ice hydrometeors (cloud ice, snow,
!           graupel, hail) to be consistent with mixing ratios using
!           diagnostic relationships.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constants (mimax, mi0, ns0, rhog, etc.) from comphy/commath.
!   - No synchronization constructs.
!   - Uses intrinsic min(), max(), sqrt() - all GPU compatible.
!   - Conditional on haiopt determines if hail is processed.
!   - All grid points are independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i) loops.
!   - Module constants can be passed as scalars to device.
! Runtime:
!   - Calls: 1080
!   - AvgLoops: 102.4M
!   - TotalTime: 13.236s (0.44%)
!   - AvgTime: 12.255ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
