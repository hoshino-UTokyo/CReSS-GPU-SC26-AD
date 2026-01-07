# Kernel 153: s_gsmoov

## Source Location
- **File**: Src/gsmoov.f90
- **Subroutine**: s_gsmoov
- **Line**: ~267
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Update y-velocity GPV data by adding diffusion term scaled by dtcoe.

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: gsmoov.f90 :: s_gsmoov
! Summary : Update y-velocity GPV data by adding diffusion term scaled by dtcoe.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global writes; updates vgpv array in-place
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
