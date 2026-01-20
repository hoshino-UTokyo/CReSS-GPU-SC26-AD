# Kernel 037: s_bcten

## Source Location
- **File**: Src/bcten.f90
- **Subroutine**: s_bcten
- **Line**: ~129

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets bottom and top boundary conditions for optional tensor array

## Runtime Profile (from test_real)
- **Calls**: 1440
- **Average Loop Length**: 811.8K
- **Total Time**: 0.068s
- **Average Time per Call**: 0.047ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bcten.f90 :: s_bcten
! Summary : Sets bottom and top boundary conditions for optional tensor array
!           by copying or negating values at boundary layers.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple array assignments to boundary planes (k=1, k=nk)
!   - No synchronization constructs beyond implicit barriers at omp end do
!   - Conditional branches based on bbc/tbc values (control flow divergence)
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Use collapse(2) for nested i,j loops to increase parallelism
! Runtime:
!   - Calls: 1440
!   - AvgLoops: 811.8K
!   - TotalTime: 0.068s (0.00%)
!   - AvgTime: 0.047ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
