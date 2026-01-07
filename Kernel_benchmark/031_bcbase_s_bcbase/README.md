# Kernel 031: s_bcbase

## Source Location
- **File**: Src/bcbase.f90
- **Subroutine**: s_bcbase
- **Line**: ~184

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Sets bottom and top boundary conditions for base state variables

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 809.1K
- **Total Time**: 0.000s
- **Average Time per Call**: 0.447ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bcbase.f90 :: s_bcbase
! Summary : Sets bottom and top boundary conditions for base state variables
!           (ubr, vbr, ptbr, qvbr, ptvbr, pibr, pbr, rbr) using extrapolation.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Multiple 2D loops over i,j for different variables
!   - Uses physical constants from m_comphy (g, cp, rd, p0)
!   - Exner function BC requires exp/log calculations
!   - Pressure and density BCs depend on previously computed pibr and ptvbr
! Next:
!   - Convert to OpenACC with collapsed i,j loops
!   - Ensure data dependencies between loops are respected (ptvbr before pibr, pibr before pbr/rbr)
!   - Consider fusing independent loops for better kernel efficiency
! Runtime:
!   - Calls: 1
!   - AvgLoops: 809.1K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.447ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
