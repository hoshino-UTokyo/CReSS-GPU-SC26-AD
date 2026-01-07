# Kernel 355: s_var8w8s

## Source Location
- **File**: Src/var8w8s.f90
- **Subroutine**: s_var8w8s
- **Line**: ~105

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Averages variable at w points to scalar points using vertical

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: var8w8s.f90 :: s_var8w8s
! Summary : Averages variable at w points to scalar points using vertical
!           2-point averaging between adjacent k levels.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global/module variable writes, only local array writes
!   - No synchronization constructs (barrier, critical, atomic)
!   - Simple loop structure with private loop indices
! Next:
!   - Direct conversion to OpenACC parallel loop or OpenACC
!   - Consider collapsing nested loops for better GPU utilization
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
