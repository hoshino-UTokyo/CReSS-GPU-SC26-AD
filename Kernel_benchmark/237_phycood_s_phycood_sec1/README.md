# Kernel 237: s_phycood

## Source Location
- **File**: Src/phycood.f90
- **Subroutine**: s_phycood
- **Line**: ~203
- **Section**: 1 of 3 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Find the highest terrain height using max reduction over 2D domain.

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 810.0K
- **Total Time**: 0.000s
- **Average Time per Call**: 0.016ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: phycood.f90 :: s_phycood
! Summary : Find the highest terrain height using max reduction over 2D domain.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses reduction(max: htmax) clause
!   - Reads from ht array, no writes to shared arrays
!   - Simple 2D loop with reduction operation
! Next:
!   - Use OpenACC parallel loop with reduction(max:htmax)
!   - GPU reductions are well supported in OpenACC
!   - May need atomic or tree-based reduction for performance
! Runtime:
!   - Calls: 1
!   - AvgLoops: 810.0K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.016ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
