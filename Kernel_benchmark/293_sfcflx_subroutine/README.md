# Kernel 293: subroutine

## Source Location
- **File**: Src/sfcflx.f90
- **Subroutine**: subroutine
- **Line**: ~161

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates exchange coefficients for surface momentum, heat,

## Runtime Profile (from test_real)
- **Calls**: 361
- **Average Loop Length**: 806.4K
- **Total Time**: 0.015s
- **Average Time per Call**: 0.043ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: sfcflx.f90 :: subroutine s_sfcflx
! Summary : Calculates exchange coefficients for surface momentum, heat,
!           and moisture fluxes from bulk coefficients.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constants (da0, dv0, wkappa) from comphy.
!   - No synchronization constructs.
!   - Simple 2D loop over (i,j) - surface arrays only.
!   - Land/sea conditional for moisture exchange coefficient.
!   - All grid points are independent (embarrassingly parallel).
!   - Uses intrinsic log() outside parallel region.
! Next:
!   - Direct OpenACC kernels should work well.
!   - 2D arrays only - good memory access pattern.
!   - Consider fusing with bulksfc call if beneficial.
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.015s (0.00%)
!   - AvgTime: 0.043ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
