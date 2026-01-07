# Kernel 046: s_buoywsi

## Source Location
- **File**: Src/buoywsi.f90
- **Subroutine**: s_buoywsi
- **Line**: ~167

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates buoyancy forcing for vertical velocity in small

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 100.4M
- **Total Time**: 143.418s
- **Average Time per Call**: 9.960ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: buoywsi.f90 :: s_buoywsi
! Summary : Calculates buoyancy forcing for vertical velocity in small
!           time step integration using implicit vertical method.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple arithmetic operations only
!   - Single conditional branch based on gwmopt
!   - Two phases: compute wb8s, then vertically average to fw
!   - Sequential dependency between phases
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Use collapse(2) for nested i,j loops
!   - Consider separate target regions for wb8s computation and fw update
!   - Simple structure well-suited for GPU offload
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 100.4M
!   - TotalTime: 143.418s (4.81%)
!   - AvgTime: 9.960ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
