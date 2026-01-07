# Kernel 126: s_getkref

## Source Location
- **File**: Src/getkref.f90
- **Subroutine**: s_getkref
- **Line**: ~183

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Determines reference index for base state pressure interpolation

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getkref.f90 :: s_getkref
! Summary : Determines reference index for base state pressure interpolation
!           by finding lowest/highest data planes that intersect model grid
! GPU diff: Hard
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic min() and max() functions
!   - Contains reduction operations (min: zdmin, altmin) (max: zdmax)
!   - Contains !$omp single regions with sequential do loops (do_kd, do_k_1, do_k_2)
!   - Module variables lim36, lim36n used from m_commath
!   - Complex control flow with early exit loops
!   - Sequential dependencies in single regions finding kdbot, kbot, ktop
! Next:
!   - Reductions can be ported but single regions need restructuring
!   - Consider splitting into separate kernels: one for reductions, one for sequential search
!   - Single regions may need to remain on host or use atomic operations
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
