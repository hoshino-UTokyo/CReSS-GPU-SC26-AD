# Kernel 387: s_xy2ll

## Source Location
- **File**: Src/xy2ll.f90
- **Subroutine**: s_xy2ll
- **Line**: ~212

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Convert x,y map coordinates to latitude/longitude using various

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 810.0K
- **Total Time**: 0.000s
- **Average Time per Call**: 0.028ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: xy2ll.f90 :: s_xy2ll
! Summary : Convert x,y map coordinates to latitude/longitude using various
!           map projections (lat-lon, Polar Stereographic, Lambert, Mercator, etc.)
! GPU diff: Medium
! Findings:
!   - No omp_get_thread usage
!   - No external function calls; uses intrinsics only (atan, cos, exp, log, sqrt)
!   - No global variable writes (only output arrays lat, lon)
!   - No explicit synchronization constructs
!   - Multiple conditional branches (mpopt) with separate omp do regions
!   - All omp do regions are mutually exclusive (only one executes per call)
!   - Element-wise computation with no loop-carried dependencies
! Next:
!   - Consider restructuring branches into separate kernels or use runtime selection
!   - Map x, y, cpj as to, and lat, lon as from
!   - Intrinsic math functions are GPU-compatible
! Runtime:
!   - Calls: 1
!   - AvgLoops: 810.0K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.028ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
