# Kernel 248: s_rbcpt

## Source Location
- **File**: Src/rbcpt.f90
- **Subroutine**: s_rbcpt
- **Line**: ~267

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Sets radiative lateral boundary conditions for potential

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rbcpt.f90 :: s_rbcpt
! Summary : Sets radiative lateral boundary conditions for potential
!           temperature perturbation at domain corners and edges (W/E/S/N).
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Multiple conditional branches based on MPI subdomain position (ebs, ebn, ebw, ebe, isub, jsub)
!   - Many small omp do regions (50+) with varying loop bounds
!   - Writes to ptpf (3D inout array) at boundary points only
!   - No sync constructs; implicit barriers at omp end do
!   - Uses module variables from m_commpi (ebs, ebn, ebw, ebe, isub, jsub, nisub, njsub)
!   - Complex conditional logic (advopt, nggopt, lspopt, vspopt options)
!   - Serial k-loop with nested parallel j or i loops in some regions
! Next:
!   - Restructure boundary updates into separate GPU kernels per edge
!   - Consider batching corner/edge updates to reduce kernel launch overhead
!   - MPI-related conditionals may require host-side decision before GPU kernel
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
