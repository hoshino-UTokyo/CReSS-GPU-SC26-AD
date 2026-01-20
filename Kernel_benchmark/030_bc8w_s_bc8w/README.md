# Kernel 030: s_bc8w

## Source Location
- **File**: Src/bc8w.f90
- **Subroutine**: s_bc8w
- **Line**: ~132

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets bottom and top boundary conditions for optional variable at w points

## Runtime Profile (from test_real)
- **Calls**: 2
- **Average Loop Length**: 811.8K
- **Total Time**: 0.000s
- **Average Time per Call**: 0.042ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bc8w.f90 :: s_bc8w
! Summary : Sets bottom and top boundary conditions for optional variable at w points
!           by copying from adjacent vertical levels based on BC type.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - 2D loops over i,j with fixed k indices (boundaries)
!   - Multiple conditional branches based on BC type (bbc, tbc)
!   - Simple array copy operations
! Next:
!   - Convert to OpenACC with collapsed i,j loops
!   - Merge bottom and top BC loops into single kernel if both are same type
! Runtime:
!   - Calls: 2
!   - AvgLoops: 811.8K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.042ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
