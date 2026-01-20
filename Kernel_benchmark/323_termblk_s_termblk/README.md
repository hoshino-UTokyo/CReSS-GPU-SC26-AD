# Kernel 323: s_termblk

## Source Location
- **File**: Src/termblk.f90
- **Subroutine**: s_termblk
- **Line**: ~305

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate terminal velocities for cloud water, rain, ice, snow,

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.4M
- **Total Time**: 13.304s
- **Average Time per Call**: 36.956ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: termblk.f90 :: s_termblk
! Summary : Calculate terminal velocities for cloud water, rain, ice, snow,
!           graupel, and hail based on microphysics options
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - Calls intrinsic functions only (exp, log, sqrt)
!   - Writes to multiple output arrays (ucq, urq, uiq, usq, ugq, uhq, ucn, urn, uin, usn, ugn, uhn)
!   - Multiple conditional branches based on flqcqi_opt, cphopt, haiopt
!   - No synchronization constructs within parallel region
! Next:
!   - Use OpenACC or OpenACC for GPU offloading
!   - Consider kernel fusion for related velocity calculations
!   - Map all input/output arrays to device
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 13.304s (0.45%)
!   - AvgTime: 36.956ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
