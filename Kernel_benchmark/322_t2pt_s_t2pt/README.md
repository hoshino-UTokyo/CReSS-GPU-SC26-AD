# Kernel 322: s_t2pt

## Source Location
- **File**: Src/t2pt.f90
- **Subroutine**: s_t2pt
- **Line**: ~118

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Convert temperature to potential temperature using pressure

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: t2pt.f90 :: s_t2pt
! Summary : Convert temperature to potential temperature using pressure
!           and Poisson equation (T * (p0/p)^(rd/cp))
! GPU diff: Easy
! Findings:
!   - Simple element-wise calculation with exp/log
!   - No conditionals within loops
!   - No function calls within parallel region
!   - Independent grid point calculations
! Next:
!   - Straightforward GPU offloading with collapse clause
!   - Consider using fast math for exp/log if accuracy permits
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
