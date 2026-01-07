# Kernel 200: s_mapfct

## Source Location
- **File**: Src/mapfct.f90
- **Subroutine**: s_mapfct
- **Line**: ~223

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate map scale factors for various projection methods (spherical,

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 810.0K
- **Total Time**: 0.001s
- **Average Time per Call**: 0.646ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: mapfct.f90 :: s_mapfct
! Summary : Calculate map scale factors for various projection methods (spherical,
!           polar stereographic, Lambert conformal, Mercator, etc.) at scalar/u/v points
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Intrinsic math functions (cos, sin, tan, exp, log, sqrt) used inside loops
!   - Multiple worksharing constructs with branching logic based on mpopt
!   - Writes to mf, mf8u, mf8v, rmf, rmf8u, rmf8v, tmp1 arrays
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Collapse nested i,j loops for better GPU occupancy
!   - Math intrinsics are GPU-compatible
! Runtime:
!   - Calls: 1
!   - AvgLoops: 810.0K
!   - TotalTime: 0.001s (0.00%)
!   - AvgTime: 0.646ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
