# Kernel 155: s_gsmoow

## Source Location
- **File**: Src/gsmoow.f90
- **Subroutine**: s_gsmoow
- **Line**: ~256
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Update z-velocity GPV data by adding diffusion term scaled by dtcoe.

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: gsmoow.f90 :: s_gsmoow
! Summary : Update z-velocity GPV data by adding diffusion term scaled by dtcoe.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global writes; updates wgpv array in-place
!   - No sync constructs
!   - Simple element-wise update operation
! Next:
!   - Convert to OpenACC with collapse(2) on j,i loops
!   - Data managed automatically via Unified Memory
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
