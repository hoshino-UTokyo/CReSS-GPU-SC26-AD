# Kernel 146: s_getzlow

## Source Location
- **File**: Src/getzlow.f90
- **Subroutine**: s_getzlow
- **Line**: ~103

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate height of lowest model level above terrain by

## Runtime Profile (from test_real)
- **Calls**: 361
- **Average Loop Length**: 806.4K
- **Total Time**: 0.012s
- **Average Time per Call**: 0.032ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getzlow.f90 :: s_getzlow
! Summary : Calculate height of lowest model level above terrain by
!           averaging vertical spacing between levels 2 and 3.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Simple 2D loop with element-wise arithmetic
!   - Reads from 3D zph array at fixed k indices (2,3)
!   - No global writes, only output array za is modified
!   - No synchronization constructs
! Next:
!   - Direct translation to OpenACC with teams distribute
!   - Consider loop collapse for j,i dimensions
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.012s (0.00%)
!   - AvgTime: 0.032ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
