# Kernel 223: s_paractl

## Source Location
- **File**: Src/paractl.f90
- **Subroutine**: s_paractl
- **Line**: ~734
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute latitude array for Mercator projection using exponential

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: paractl.f90 :: s_paractl
! Summary : Compute latitude array for Mercator projection using exponential
!           and trigonometric functions.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic functions: atan, exp, max, min (GPU-compatible)
!   - No global/module variable writes (only mlat output array)
!   - No sync constructs
!   - Conditional based on uniopt_uni with different loop bounds
! Next:
!   - Convert to OpenACC with device math library
!   - Ensure atan/exp are available on GPU device
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
