# Kernel 386: s_xy2ij

## Source Location
- **File**: Src/xy2ij.f90
- **Subroutine**: s_xy2ij
- **Line**: ~150

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Convert 2D x,y coordinates to real grid indices (ri, rj) using

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: xy2ij.f90 :: s_xy2ij
! Summary : Convert 2D x,y coordinates to real grid indices (ri, rj) using
!           inverse grid spacing multipliers.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - No global variable writes (only output arrays ri, rj)
!   - No synchronization constructs beyond implicit barrier at end do
!   - Simple element-wise computation with no data dependencies
! Next:
!   - Direct translation to OpenACC or OpenACC parallel loop
!   - Map x2d, y2d as to, and ri, rj as from
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
