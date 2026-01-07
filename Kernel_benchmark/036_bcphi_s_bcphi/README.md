# Kernel 036: s_bcphi

## Source Location
- **File**: Src/bcphi.f90
- **Subroutine**: s_bcphi
- **Line**: ~114

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets bottom and top boundary conditions for parabolic PDE solver

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bcphi.f90 :: s_bcphi
! Summary : Sets bottom and top boundary conditions for parabolic PDE solver
!           by copying from adjacent vertical levels.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Simple 2D loop over i,j with fixed k indices
!   - No external module dependencies inside parallel region
!   - Straightforward array copy operations
! Next:
!   - Convert to OpenACC with collapsed i,j loops
!   - Very simple kernel, good candidate for GPU porting
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
