# Kernel 142: s_getvdens

## Source Location
- **File**: Src/getvdens.f90
- **Subroutine**: s_getvdens
- **Line**: ~104

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute inverse of base state density (1/rbr) for all 3D

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.4M
- **Total Time**: 1.102s
- **Average Time per Call**: 3.061ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getvdens.f90 :: s_getvdens
! Summary : Compute inverse of base state density (1/rbr) for all 3D
!           grid points for use in momentum equations.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Simple element-wise division: rbv = 1/rbr
!   - No global writes, only output array rbv is modified
!   - No synchronization constructs other than implicit barrier
! Next:
!   - Direct translation to OpenACC with collapsed loops
!   - Consider loop collapse for k,j,i dimensions
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 1.102s (0.04%)
!   - AvgTime: 3.061ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
