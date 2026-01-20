# Kernel 171: s_hintrdr

## Source Location
- **File**: Src/hintrdr.f90
- **Subroutine**: s_hintrdr
- **Line**: ~184

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate grid distances and perform horizontal interpolation

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: hintrdr.f90 :: s_hintrdr
! Summary : Calculate grid distances and perform horizontal interpolation
!           for radar data with missing value handling (lim34n, lim35n).
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region; uses intrinsics only
!   - Distance calculation for fproc='cal', then interpolation
!   - Complex branching for missing value checks (lim34n/lim35n)
!   - Different code paths for mpopt < 10 vs >= 10
!   - Outer k-loop with inner j,i loops parallelized
!   - Writes to di, dj (distances) and var (interpolated values)
!   - No sync constructs
! Next:
!   - Convert to OpenACC with collapse(2) or collapse(3)
!   - Missing value conditionals may cause GPU thread divergence
!   - Consider masking approach for missing values
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
