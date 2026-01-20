# Kernel 207: s_newblk_noevap

## Source Location
- **File**: Src/newblk_noevap.f90
- **Subroutine**: s_newblk_noevap
- **Line**: ~369

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Solve new potential temperature perturbation, mixing ratios, and

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: newblk_noevap.f90 :: s_newblk_noevap
! Summary : Solve new potential temperature perturbation, mixing ratios, and
!           concentrations for bulk microphysics without temperature variation for evaporation
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic abs() and max() functions - GPU compatible
!   - Uses module constants from m_comphy (cp, mr0, mi0, ms0)
!   - Very complex conditional logic based on cphopt (2, 3, or 4)
!   - Different calculations for nk==1 (2D) vs nk>1 (3D)
!   - Many private variables with complex local computations
!   - Writes to ptpf, qvf, qcf, qrf, qif, qsf, qgf, nccf, ncrf, ncif, ncsf, ncgf
!   - Multiple nested conditionals checking thresq thresholds
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Consider restructuring to reduce branch divergence on GPU
!   - May benefit from separating cphopt cases into different kernels
!   - Collapse nested i,j loops for better GPU occupancy
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
