# Kernel 137: s_getta2d

## Source Location
- **File**: Src/getta2d.f90
- **Subroutine**: s_getta2d
- **Line**: ~106

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates 2D air temperature from potential temperature and

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getta2d.f90 :: s_getta2d
! Summary : Calculates 2D air temperature from potential temperature and
!           Exner function: t = (ptbr + ptp) * pi
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Simple element-wise multiplication at each grid point
!   - 2D loop over surface grid points (i,j)
!   - No loop-carried dependencies
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Direct port to OpenACC parallel loop
!   - Collapse j,i loops for better GPU occupancy
!   - Trivial computation, ensure data is already on device
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
