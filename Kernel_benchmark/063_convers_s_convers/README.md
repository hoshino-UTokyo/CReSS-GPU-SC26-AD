# Kernel 063: s_convers

## Source Location
- **File**: Src/convers.f90
- **Subroutine**: s_convers
- **Line**: ~228

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates conversion rates between cloud/ice phases: cloud water

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 4.885s
- **Average Time per Call**: 0.107ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: convers.f90 :: s_convers
! Summary : Calculates conversion rates between cloud/ice phases: cloud water
!           to rain, cloud ice to snow, and snow to graupel.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls (only intrinsic abs, exp, log, min, sqrt)
!   - Complex conditional logic with cphopt branching (abs(cphopt)==2 vs >=3)
!   - Uses module variables from m_commath and m_comphy (read-only constants)
!   - Separate paths for nk=1 and nk>1 cases
!   - Multiple output arrays (cncr, cnis, cnsg, cnsgn)
! Next:
!   - Can be ported with OpenACC with loop collapse
!   - May benefit from separating cphopt==2 and cphopt>=3 into distinct kernels
!   - Consider constant memory for module physical constants
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 4.885s (0.16%)
!   - AvgTime: 0.107ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
