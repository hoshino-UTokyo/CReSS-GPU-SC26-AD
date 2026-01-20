# Kernel 044: s_buoywb

## Source Location
- **File**: Src/buoywb.f90
- **Subroutine**: s_buoywb
- **Line**: ~164

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates buoyancy forcing for vertical velocity in large

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 100.4M
- **Total Time**: 5.412s
- **Average Time per Call**: 15.032ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: buoywb.f90 :: s_buoywb
! Summary : Calculates buoyancy forcing for vertical velocity in large
!           time step integration for dry/moist air with/without microphysics.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Multiple conditional branches based on fmois, gwmopt, cphopt
!   - Three phases: compute qvd, compute wb8s, vertically average to wfrc
!   - Simple arithmetic but significant control flow divergence
!   - Sequential dependencies between phases (qvd -> wb8s -> wfrc)
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Use collapse(2) for nested i,j loops
!   - Consider separate target regions for each major conditional branch
!   - Ensure proper data movement for qvd, wb8s work arrays
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.4M
!   - TotalTime: 5.412s (0.18%)
!   - AvgTime: 15.032ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
