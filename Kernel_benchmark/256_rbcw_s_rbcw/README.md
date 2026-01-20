# Kernel 256: s_rbcw

## Source Location
- **File**: Src/rbcw.f90
- **Subroutine**: s_rbcw
- **Line**: ~247

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Sets radiative lateral boundary conditions for z-velocity (w)

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 126
- **Total Time**: 1.864s
- **Average Time per Call**: 0.129ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rbcw.f90 :: s_rbcw
! Summary : Sets radiative lateral boundary conditions for z-velocity (w)
!           at domain corners and edges (W/E/S/N) with phase speed updates.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Multiple conditional branches based on MPI subdomain position (ebs, ebn, ebw, ebe, isub, jsub)
!   - Many omp do regions: 8 for corners + 8 for edges = 16 total
!   - Writes to w (3D inout array) at boundary points, in-place update
!   - No sync constructs; implicit barriers at omp end do
!   - Uses module variables from m_commpi for domain decomposition
!   - Conditional on gpvvar, nggopt, lspopt, vspopt for GPV nudging
!   - Corner updates involve combined x and y phase speeds (radwe, radsn)
! Next:
!   - Separate corner and edge kernels for GPU
!   - MPI conditionals evaluated on host before kernel launch
!   - Consider batching corner updates to reduce kernel overhead
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 126
!   - TotalTime: 1.864s (0.06%)
!   - AvgTime: 0.129ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
