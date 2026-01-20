# Kernel 152: s_gsmoov

## Source Location
- **File**: Src/gsmoov.f90
- **Subroutine**: s_gsmoov
- **Line**: ~159
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute 2nd-order diffusion term for y-velocity (vgpv) smoothing

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: gsmoov.f90 :: s_gsmoov
! Summary : Compute 2nd-order diffusion term for y-velocity (vgpv) smoothing
!           using 3D stencil in x, y, z directions.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global writes; only writes to dfv array (output buffer)
!   - No sync constructs (barriers, critical sections)
!   - Simple stencil computation with k-loop parallelized
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
