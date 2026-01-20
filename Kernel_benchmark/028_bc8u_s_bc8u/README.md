# Kernel 028: s_bc8u

## Source Location
- **File**: Src/bc8u.f90
- **Subroutine**: s_bc8u
- **Line**: ~134

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets west and east boundary conditions for optional variable at u points

## Runtime Profile (from test_real)
- **Calls**: 720
- **Average Loop Length**: 115.3K
- **Total Time**: 1.712s
- **Average Time per Call**: 2.378ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bc8u.f90 :: s_bc8u
! Summary : Sets west and east boundary conditions for optional variable at u points
!           by copying from adjacent interior points based on BC type.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Nested loops with outer k-loop serial, inner j-loop parallelized
!   - Uses module variables from m_commpi (ebw, ebe, isub, nisub)
!   - Conditional execution based on BC type (wbc, ebc) and subdomain position
! Next:
!   - Convert to OpenACC with collapsed j,k loops
!   - Restructure loops to have k as inner loop for better GPU coalescing
! Runtime:
!   - Calls: 720
!   - AvgLoops: 115.3K
!   - TotalTime: 1.712s (0.06%)
!   - AvgTime: 2.378ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
