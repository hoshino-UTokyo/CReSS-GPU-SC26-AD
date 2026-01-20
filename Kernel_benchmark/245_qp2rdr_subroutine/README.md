# Kernel 245: subroutine

## Source Location
- **File**: Src/qp2rdr.f90
- **Subroutine**: subroutine
- **Line**: ~144

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Performs analysis nudging of precipitation mixing ratio

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: qp2rdr.f90 :: subroutine s_qp2rdr
! Summary : Performs analysis nudging of precipitation mixing ratio
!           toward radar observations with time interpolation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constant lim34n from commath.
!   - No synchronization constructs.
!   - Uses intrinsic max() - GPU compatible.
!   - Conditional updates based on data validity thresholds.
!   - All grid points are independent.
! Next:
!   - Direct OpenACC kernels with conditional update.
!   - Data validity check may cause minor warp divergence.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
