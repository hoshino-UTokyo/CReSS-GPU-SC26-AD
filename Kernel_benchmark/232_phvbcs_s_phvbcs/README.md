# Kernel 232: s_phvbcs

## Source Location
- **File**: Src/phvbcs.f90
- **Subroutine**: s_phvbcs
- **Line**: ~302

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Calculate differential phase speed terms for scalar variables

## Runtime Profile (from test_real)
- **Calls**: 1080
- **Average Loop Length**: 112.2K
- **Total Time**: 11.834s
- **Average Time per Call**: 10.958ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: phvbcs.f90 :: s_phvbcs
! Summary : Calculate differential phase speed terms for scalar variables
!           at west/east/south/north boundaries for radiation BCs.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region (uses intrinsics: abs, max, min, sign)
!   - No global variable writes (outputs to scpx, scpy boundary arrays)
!   - No explicit sync constructs
!   - Complex conditional logic based on boundary condition types (wbc, ebc, sbc, nbc)
!   - Uses MPI domain decomposition variables (ebw, ebe, ebs, ebn, isub, jsub)
!   - Reduction-like pattern for cpavex, cpavey with nkm3v scaling
!   - Many conditional branches affecting control flow
! Next:
!   - Consider GPU porting only for large domains where boundary computation is significant
!   - Multiple kernel launches may be needed for different boundary conditions
!   - Reduction operations for cpavex, cpavey need GPU reduction support
!   - MPI rank checks (isub, jsub) determine which boundaries are active
! Runtime:
!   - Calls: 1080
!   - AvgLoops: 112.2K
!   - TotalTime: 11.834s (0.40%)
!   - AvgTime: 10.958ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
