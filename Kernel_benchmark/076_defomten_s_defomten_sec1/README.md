# Kernel 076: s_defomten

## Source Location
- **File**: Src/defomten.f90
- **Subroutine**: s_defomten
- **Line**: ~277
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Initialize velocity arrays at w-points for terrain-following

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: defomten.f90 :: s_defomten (first parallel region)
! Summary : Initialize velocity arrays at w-points for terrain-following
!           coordinate deformation tensor calculation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Private variable k for outer loop
!   - Writes to s11, s22 arrays (temporary storage for u, v at w-points)
!   - Simple stencil averaging in vertical direction
! Next:
!   - Direct conversion to OpenACC with data region
!   - Can collapse k and j loops for better GPU occupancy
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
