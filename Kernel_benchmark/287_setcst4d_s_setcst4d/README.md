# Kernel 287: s_setcst4d

## Source Location
- **File**: Src/setcst4d.f90
- **Subroutine**: s_setcst4d
- **Line**: ~126

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Fill 4D array with a constant value (real type)

## Runtime Profile (from test_real)
- **Calls**: 25
- **Average Loop Length**: 158.7M
- **Total Time**: 0.797s
- **Average Time per Call**: 31.877ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: setcst4d.f90 :: s_setcst4d
! Summary : Fill 4D array with a constant value (real type)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls
!   - Simple 4D loop with direct assignment
!   - No synchronization constructs
!   - Trivially parallelizable
! Next:
!   - Convert to OpenACC with parallel loop collapse(4)
!   - Consider using memset or array assignment for better performance
! Runtime:
!   - Calls: 25
!   - AvgLoops: 158.7M
!   - TotalTime: 0.797s (0.03%)
!   - AvgTime: 31.877ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
