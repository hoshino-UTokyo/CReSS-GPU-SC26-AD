# Kernel 247: s_radiat

## Source Location
- **File**: Src/radiat.f90
- **Subroutine**: s_radiat
- **Line**: ~300

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates zenith angle, short/long wave radiation fluxes

## Runtime Profile (from test_real)
- **Calls**: 361
- **Average Loop Length**: 806.4K
- **Total Time**: 1.831s
- **Average Time per Call**: 5.071ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: radiat.f90 :: s_radiat
! Summary : Calculates zenith angle, short/long wave radiation fluxes
!           (rgd, rsd, rld, rlu) based on dry/moist air conditions.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - Intrinsic functions used: cos, sin, exp, log10, sqrt, max, min
!   - Multiple omp do regions inside single parallel region
!   - Writes to zph8s, zref, coseta (2D), rgd, rsd, rld, rlu (2D output arrays)
!   - No sync constructs; implicit barriers at omp end do
!   - Conditional branching based on fmois (dry/moist) and cphopt
!   - Serial k-loop with nested parallel i,j loops
!   - Uses module constants from m_comdays, m_commath, m_comphy
! Next:
!   - Collapse k-loop with i,j loops if possible
!   - Consider separating dry/moist code paths for GPU kernels
!   - Hoist conditional checks outside parallel region if feasible
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 1.831s (0.06%)
!   - AvgTime: 5.071ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
