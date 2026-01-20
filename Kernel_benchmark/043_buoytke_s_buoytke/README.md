# Kernel 043: s_buoytke

## Source Location
- **File**: Src/buoytke.f90
- **Subroutine**: s_buoytke
- **Line**: ~128

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates buoyancy production term for turbulent kinetic

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 101.2M
- **Total Time**: 3.027s
- **Average Time per Call**: 8.409ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: buoytke.f90 :: s_buoytke
! Summary : Calculates buoyancy production term for turbulent kinetic
!           energy equation using Brunt-Vaisala frequency and diffusivity.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple arithmetic operations only
!   - Two sequential loop nests: first computes tmp1, second updates tkefrc
!   - Data dependency: tmp1 must be computed before tkefrc update
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Use collapse(2) for nested i,j loops within each k-loop
!   - Keep two separate target regions or use explicit barrier between phases
!   - Consider fusing loops if tmp1 dependency can be restructured
! Runtime:
!   - Calls: 360
!   - AvgLoops: 101.2M
!   - TotalTime: 3.027s (0.10%)
!   - AvgTime: 8.409ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
