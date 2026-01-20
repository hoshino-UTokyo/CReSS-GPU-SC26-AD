# Kernel 371: s_vint133g

## Source Location
- **File**: Src/vint133g.f90
- **Subroutine**: s_vint133g
- **Line**: ~161

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Interpolate horizontally-varying variable to model grid (generic version)

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vint133g.f90 :: s_vint133g
! Summary : Interpolate horizontally-varying variable to model grid (generic version)
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - getindx() called before parallel region (safe)
!   - Reads from zph, invar (3D), z1d (1D); writes to outvar (3D)
!   - First section: extrapolation with k loop serial, i,j parallelized
!   - Second section: kl,k loops serial, i,j parallelized for interpolation
!   - Variable arrangement flag (xo) determines loop bounds
! Next:
!   - Collapse k,j,i loops for GPU parallelism
!   - Use OpenACC teams distribute parallel do collapse(3)
!   - Consider restructuring to compute kl index per (i,j,k) point directly
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
