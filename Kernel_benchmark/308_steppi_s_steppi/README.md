# Kernel 308: s_steppi

## Source Location
- **File**: Src/steppi.f90
- **Subroutine**: s_steppi
- **Line**: ~215

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Advances pressure perturbation in time using forcing term

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 100.4M
- **Total Time**: 62.933s
- **Average Time per Call**: 4.370ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: steppi.f90 :: s_steppi
! Summary : Advances pressure perturbation in time using forcing term
!           with horizontally explicit/vertically implicit method
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Single !$omp do loop with schedule(runtime)
!   - Simple element-wise update: ppf += dts*fp/jcb
!   - Writes only to ppf array
!   - No synchronization constructs besides implicit barriers
! Next:
!   - Data managed automatically via Unified Memory
!   - Convert to !$acc parallel loop collapse(2)
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 100.4M
!   - TotalTime: 62.933s (2.11%)
!   - AvgTime: 4.370ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
