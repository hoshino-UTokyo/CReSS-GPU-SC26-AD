# Kernel 075: s_defomssq

## Source Location
- **File**: Src/defomssq.f90
- **Subroutine**: s_defomssq
- **Line**: ~125

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate magnitude of deformation tensor squared from diagonal

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.4M
- **Total Time**: 3.000s
- **Average Time per Call**: 8.333ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: defomssq.f90 :: s_defomssq
! Summary : Calculate magnitude of deformation tensor squared from diagonal
!           and off-diagonal strain rate components using stencil averaging.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Private variable k for outer loop
!   - Writes to ssq output array
!   - Stencil operations averaging s12, s31, s32 at neighboring points
!   - Independent operations for each grid point
! Next:
!   - Direct conversion to OpenACC with collapsed loops
!   - Data managed automatically via Unified Memory
!   - Good candidate for GPU due to arithmetic intensity
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 3.000s (0.10%)
!   - AvgTime: 8.333ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
