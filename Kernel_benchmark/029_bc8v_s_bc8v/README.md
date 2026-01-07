# Kernel 029: s_bc8v

## Source Location
- **File**: Src/bc8v.f90
- **Subroutine**: s_bc8v
- **Line**: ~134

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets south and north boundary conditions for optional variable at v points

## Runtime Profile (from test_real)
- **Calls**: 720
- **Average Loop Length**: 115.3K
- **Total Time**: 1.731s
- **Average Time per Call**: 2.404ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bc8v.f90 :: s_bc8v
! Summary : Sets south and north boundary conditions for optional variable at v points
!           by copying from adjacent interior points based on BC type.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Nested loops with outer k-loop serial, inner i-loop parallelized
!   - Uses module variables from m_commpi (ebs, ebn, jsub, njsub)
!   - Conditional execution based on BC type (sbc, nbc) and subdomain position
! Next:
!   - Convert to OpenACC with collapsed i,k loops
!   - Restructure loops to have k as inner loop for better GPU coalescing
! Runtime:
!   - Calls: 720
!   - AvgLoops: 115.3K
!   - TotalTime: 1.731s (0.06%)
!   - AvgTime: 2.404ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
