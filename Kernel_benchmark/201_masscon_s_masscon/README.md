# Kernel 201: s_masscon

## Source Location
- **File**: Src/masscon.f90
- **Subroutine**: s_masscon
- **Line**: ~300

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate the known quantity (RHS) for mass consistent velocity fitting

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: masscon.f90 :: s_masscon
! Summary : Calculate the known quantity (RHS) for mass consistent velocity fitting
!           by computing divergence of momentum from Jacobian-weighted velocities
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to tmp1, tmp2, tmp3, known arrays
!   - Multiple worksharing constructs with different k-loop ranges
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Collapse nested i,j loops for better GPU occupancy
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
