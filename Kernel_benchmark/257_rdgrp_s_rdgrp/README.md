# Kernel 257: s_rdgrp

## Source Location
- **File**: Src/rdgrp.f90
- **Subroutine**: s_rdgrp
- **Line**: ~451

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Find min/max indices of active group domains to determine reductional domain bounds

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 1
- **Total Time**: 0.000s
- **Average Time per Call**: 0.014ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rdgrp.f90 :: s_rdgrp
! Summary : Find min/max indices of active group domains to determine reductional domain bounds
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses min/max reductions on iwred, jsred, iered, jnred
!   - No synchronization constructs
!   - Simple 2D loop with conditional bounds checking
! Next:
!   - Convert to OpenACC with teams distribute parallel for and reduction clause
!   - Alternatively use OpenACC with parallel loop reduction
! Runtime:
!   - Calls: 1
!   - AvgLoops: 1
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.014ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
