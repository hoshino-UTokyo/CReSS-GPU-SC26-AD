# Kernel 070: s_coriuvw

## Source Location
- **File**: Src/coriuvw.f90
- **Subroutine**: s_coriuvw
- **Line**: ~136

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates Coriolis force contributions to u, v, w velocity

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: coriuvw.f90 :: s_coriuvw
! Summary : Calculates Coriolis force contributions to u, v, w velocity
!           forcing terms using staggered grid interpolations.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls
!   - Multiple separate loop nests for different components (ufrc, vfrc, wfrc)
!   - Uses temporary arrays tmp1, tmp2 for intermediate calculations
!   - Data dependency: tmp1 computed first, then used for tmp2 and forces
!   - Multiple implicit barriers between loop nests
! Next:
!   - Consider fusing loops where possible to reduce kernel launches
!   - Temporary arrays tmp1, tmp2 need to be on GPU
!   - May need explicit synchronization between kernel sections
!   - Collapse k,j,i loops within each section for better occupancy
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
