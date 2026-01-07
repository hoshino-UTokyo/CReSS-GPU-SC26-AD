# Kernel 083: s_diagni

## Source Location
- **File**: Src/diagni.f90
- **Subroutine**: s_diagni
- **Line**: ~174

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate diagnostic concentrations for all ice hydrometeor

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 102.4M
- **Total Time**: 0.012s
- **Average Time per Call**: 11.720ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diagni.f90 :: s_diagni
! Summary : Calculate diagnostic concentrations for all ice hydrometeor
!           categories (cloud ice, snow, graupel, hail) from mixing ratios.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic max, min, sqrt functions (GPU-compatible)
!   - Private variable k for outer loop; rbv for local scalar
!   - Writes to nidia output array for multiple ice categories
!   - Conditional branch based on haiopt (3 vs 4 categories)
!   - Independent point-wise operations per grid cell
! Next:
!   - Direct conversion to OpenACC with collapsed loops
!   - Handle haiopt conditional outside kernel or use single kernel with masking
!   - Data managed automatically via Unified Memory
! Runtime:
!   - Calls: 1
!   - AvgLoops: 102.4M
!   - TotalTime: 0.012s (0.00%)
!   - AvgTime: 11.720ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
