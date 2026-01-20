# Kernel 084: s_diagnsg

## Source Location
- **File**: Src/diagnsg.f90
- **Subroutine**: s_diagnsg
- **Line**: ~171

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate diagnostic concentrations for precipitation ice

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diagnsg.f90 :: s_diagnsg
! Summary : Calculate diagnostic concentrations for precipitation ice
!           categories (snow, graupel, hail) from mixing ratios and density.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic max, min, sqrt functions (GPU-compatible)
!   - Private variable k for outer loop
!   - Writes to nidia output array for snow, graupel, hail categories
!   - Conditional branch based on haiopt (2 vs 3 precipitation categories)
!   - Independent point-wise operations per grid cell
! Next:
!   - Direct conversion to OpenACC with collapsed loops
!   - Handle haiopt conditional outside kernel or use unified kernel
!   - Data managed automatically via Unified Memory
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
