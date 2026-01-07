# Kernel 365: s_vculs

## Source Location
- **File**: Src/vculs.f90
- **Subroutine**: s_vculs
- **Line**: ~171

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Compute vertical scalar advection using Cubic Lagrange scheme

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vculs.f90 :: s_vculs
! Summary : Compute vertical scalar advection using Cubic Lagrange scheme
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Reads from sp, wc8s; writes to sf (all 3D arrays)
!   - Outer k loop is serial with nested !$omp do for i,j
!   - Conditional branches (wc8s > 0) for upwind/downwind stencil selection
!   - Private variables: k (shared among worksharing), i,j,a,b,c (private per iteration)
! Next:
!   - Collapse k,j,i loops for GPU parallelism
!   - Use OpenACC teams distribute parallel do collapse(3)
!   - Data managed automatically via Unified Memory
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
