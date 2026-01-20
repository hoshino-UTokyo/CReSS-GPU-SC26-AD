# Kernel 314: s_stepwe

## Source Location
- **File**: Src/stepwe.f90
- **Subroutine**: s_stepwe
- **Line**: ~259

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Update z-velocity component using forcing and acoustic terms

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: stepwe.f90 :: s_stepwe
! Summary : Update z-velocity component using forcing and acoustic terms
!           with explicit time integration
! GPU diff: Easy
! Findings:
!   - Simple array update with element-wise operations
!   - No function calls within parallel region
!   - No conditionals within inner loops
!   - Division by rst8w for density weighting
! Next:
!   - Straightforward GPU offloading with collapse clause
!   - Consider fusing with boundary condition operations
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
