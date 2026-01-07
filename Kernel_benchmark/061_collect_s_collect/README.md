# Kernel 061: s_collect

## Source Location
- **File**: Src/collect.f90
- **Subroutine**: s_collect
- **Line**: ~333

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates collection rates between various hydrometeor species

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 898
- **Total Time**: 12.919s
- **Average Time per Call**: 0.283ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: collect.f90 :: s_collect
! Summary : Calculates collection rates between various hydrometeor species
!           (cloud water, rain, ice, snow, graupel) for bulk microphysics.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls (only intrinsic abs, exp, log, sqrt)
!   - Complex conditional logic with many if-else branches based on cphopt
!   - Multiple output arrays written independently per grid point
!   - Uses module variables from m_commath and m_comphy (read-only constants)
!   - Large number of private variables for each work item
!   - nk=1 case handled separately from nk>1 case
! Next:
!   - Consider separate kernels for cphopt==2 and cphopt>=3 cases
!   - May need register pressure optimization due to many local variables
!   - Collapse loops for better GPU utilization
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 898
!   - TotalTime: 12.919s (0.43%)
!   - AvgTime: 0.283ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
