# Kernel 105: subroutine

## Source Location
- **File**: Src/exbcu.f90
- **Subroutine**: subroutine
- **Line**: ~213

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Forces lateral boundary values of u velocity to external

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 112.2K
- **Total Time**: 3.626s
- **Average Time per Call**: 0.252ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: exbcu.f90 :: subroutine s_exbcu
! Summary : Forces lateral boundary values of u velocity to external
!           (GPV) boundary values with radiative/relaxation approach.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module variables (ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub).
!   - No synchronization constructs.
!   - Uses intrinsic abs() - GPU compatible.
!   - Boundary-position-dependent conditionals (W, E, S, N edges).
!   - Only boundary cells are updated - sparse computation.
! Next:
!   - Consider separate kernels for each boundary region.
!   - Boundary-only work has low arithmetic intensity on GPU.
!   - Domain decomposition flags need proper handling on GPU.
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 112.2K
!   - TotalTime: 3.626s (0.12%)
!   - AvgTime: 0.252ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
