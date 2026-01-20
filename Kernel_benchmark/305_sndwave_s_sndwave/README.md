# Kernel 305: s_sndwave

## Source Location
- **File**: Src/sndwave.f90
- **Subroutine**: s_sndwave
- **Line**: ~117

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Computes base state density times sound wave speed squared

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 102.4M
- **Total Time**: 0.003s
- **Average Time per Call**: 3.417ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: sndwave.f90 :: s_sndwave
! Summary : Computes base state density times sound wave speed squared
!           from base state pressure (rcsq = cp/cv * pbr)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Single !$omp do loop with schedule(runtime)
!   - Simple element-wise computation with scalar cpdvcv
!   - Writes only to rcsq array
!   - No synchronization constructs besides implicit barriers
! Next:
!   - Data managed automatically via Unified Memory
!   - Convert to !$acc parallel loop collapse(3)
! Runtime:
!   - Calls: 1
!   - AvgLoops: 102.4M
!   - TotalTime: 0.003s (0.00%)
!   - AvgTime: 3.417ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
