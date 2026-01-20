# Kernel 020: s_allociot

## Source Location
- **File**: Src/allociot.f90
- **Subroutine**: s_allociot
- **Line**: ~223

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Initialize the I/O unit number table (iolst) with sequential

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 2
- **Total Time**: 0.000s
- **Average Time per Call**: 0.061ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: allociot.f90 :: s_allociot
! Summary : Initialize the I/O unit number table (iolst) with sequential
!           unit numbers starting from 11
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to module-level array iolst from m_comionum
!   - Simple 1D loop with no data dependencies
!   - Very small loop size (nio typically small)
! Next:
!   - Trivial GPU port but likely not worth offloading
!   - Small array size means CPU execution is faster
!   - Keep on CPU; initialization is one-time cost
! Runtime:
!   - Calls: 1
!   - AvgLoops: 2
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.061ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
