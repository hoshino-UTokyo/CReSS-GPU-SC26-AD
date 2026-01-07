# Kernel 125: s_getexner

## Source Location
- **File**: Src/getexner.f90
- **Subroutine**: s_getexner
- **Line**: ~127

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates total pressure (p = pbr + pp) and Exner function

## Runtime Profile (from test_real)
- **Calls**: 1080
- **Average Loop Length**: 102.4M
- **Total Time**: 6.873s
- **Average Time per Call**: 6.364ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getexner.f90 :: s_getexner
! Summary : Calculates total pressure (p = pbr + pp) and Exner function
!           (pi = (p/p0)^(rd/cp)) at each grid point
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic exp() and log() functions - GPU compatible
!   - Simple element-wise computation, no dependencies between iterations
!   - Module constants rd, cp, p0 used from m_comphy
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Direct port to OpenACC parallel loop
!   - Ensure module constants are accessible on device
!   - Consider collapsing k,j,i loops for better GPU occupancy
! Runtime:
!   - Calls: 1080
!   - AvgLoops: 102.4M
!   - TotalTime: 6.873s (0.23%)
!   - AvgTime: 6.364ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
