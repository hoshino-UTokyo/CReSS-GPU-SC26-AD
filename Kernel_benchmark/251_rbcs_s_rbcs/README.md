# Kernel 251: s_rbcs

## Source Location
- **File**: Src/rbcs.f90
- **Subroutine**: s_rbcs
- **Line**: ~224

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Sets radiative lateral boundary conditions for optional scalar

## Runtime Profile (from test_real)
- **Calls**: 1080
- **Average Loop Length**: 125
- **Total Time**: 0.202s
- **Average Time per Call**: 0.187ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rbcs.f90 :: s_rbcs
! Summary : Sets radiative lateral boundary conditions for optional scalar
!           variable at domain corners and edges (no non-negative clamping).
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Multiple conditional branches based on MPI subdomain position (ebs, ebn, ebw, ebe, isub, jsub)
!   - Multiple omp do regions with k or (j,k)/(i,k) loop nests
!   - Writes to sf (3D inout array) at boundary points only
!   - No max() clamping unlike rbcq/rbcqv (allows negative values)
!   - No sync constructs; implicit barriers at omp end do
!   - Uses module variables from m_commpi for domain decomposition
!   - Simpler conditional logic than rbcq (no gpvvar check, no nggopt)
! Next:
!   - Separate boundary kernels for GPU (one per edge/corner)
!   - MPI conditionals should be evaluated on host before kernel launch
!   - Can share kernel structure with rbcs0 but without max() clamping
! Runtime:
!   - Calls: 1080
!   - AvgLoops: 125
!   - TotalTime: 0.202s (0.01%)
!   - AvgTime: 0.187ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
