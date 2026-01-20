# Kernel 331: s_totalqwi

## Source Location
- **File**: Src/totalqwi.f90
- **Subroutine**: s_totalqwi
- **Line**: ~135

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate total water and ice mixing ratio by summing

## Runtime Profile (from test_real)
- **Calls**: 720
- **Average Loop Length**: 102.4M
- **Total Time**: 4.982s
- **Average Time per Call**: 6.920ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: totalqwi.f90 :: s_totalqwi
! Summary : Calculate total water and ice mixing ratio by summing
!           water and ice hydrometeor categories based on cphopt/haiopt
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls
!   - Single output array (qall)
!   - Multiple conditional branches based on cphopt, haiopt
!   - Loop over bin categories for cphopt >= 11
!   - No synchronization constructs
! Next:
!   - Straightforward GPU port with conditional branches
!   - Consider specialized kernels for bulk vs bin microphysics
!   - Bin category loops can be unrolled or parallelized
! Runtime:
!   - Calls: 720
!   - AvgLoops: 102.4M
!   - TotalTime: 4.982s (0.17%)
!   - AvgTime: 6.920ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
