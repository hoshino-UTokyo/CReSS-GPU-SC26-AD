# Kernel 032: s_bcgsms

## Source Location
- **File**: Src/bcgsms.f90
- **Subroutine**: s_bcgsms
- **Line**: ~152

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Sets boundary conditions for scalar diffusion term in GPV smoothing

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bcgsms.f90 :: s_bcgsms
! Summary : Sets boundary conditions for scalar diffusion term in GPV smoothing
!           at west, east, south, north, bottom, and top boundaries.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Nested loops with outer k-loop serial, inner i or j loop parallelized
!   - Uses module variables from m_commpi (ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub)
!   - Conditional execution based on subdomain position and BC type
!   - Bottom/top BC uses 2D loop over i,j
! Next:
!   - Convert to OpenACC with collapsed loops
!   - Restructure nested loops to expose more parallelism in k dimension
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
