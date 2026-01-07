# Kernel 249: s_rbcq

## Source Location
- **File**: Src/rbcq.f90
- **Subroutine**: s_rbcq
- **Line**: ~281

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Sets radiative lateral boundary conditions for optional mixing

## Runtime Profile (from test_real)
- **Calls**: 1800
- **Average Loop Length**: 125
- **Total Time**: 0.338s
- **Average Time per Call**: 0.188ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rbcq.f90 :: s_rbcq
! Summary : Sets radiative lateral boundary conditions for optional mixing
!           ratio at domain corners and edges with non-negative clamping.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Multiple conditional branches based on MPI subdomain position (ebs, ebn, ebw, ebe, isub, jsub)
!   - Many omp do regions with k or (j,k)/(i,k) loop nests
!   - Uses max() intrinsic to clamp values >= 0
!   - Writes to qf (3D inout array) at boundary points only
!   - No sync constructs; implicit barriers at omp end do
!   - Uses module variables from m_commpi for domain decomposition
!   - Complex conditional logic based on gpvvar, advopt, nggopt, lspopt, vspopt
! Next:
!   - Separate boundary kernels for GPU (one per edge/corner)
!   - MPI conditionals should be evaluated on host before kernel launch
!   - Consider using atomic operations if boundaries overlap in GPU version
! Runtime:
!   - Calls: 1800
!   - AvgLoops: 125
!   - TotalTime: 0.338s (0.01%)
!   - AvgTime: 0.188ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
