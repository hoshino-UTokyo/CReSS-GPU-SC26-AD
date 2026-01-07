# Kernel 115: subroutine

## Source Location
- **File**: Src/gaussel.f90
- **Subroutine**: subroutine
- **Line**: ~171

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Solves tridiagonal linear systems using Gauss elimination

## Runtime Profile (from test_real)
- **Calls**: 14418
- **Average Loop Length**: 802.8K
- **Total Time**: 136.847s
- **Average Time per Call**: 9.491ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: gaussel.f90 :: subroutine s_gaussel
! Summary : Solves tridiagonal linear systems using Gauss elimination
!           (Thomas algorithm) or partial pivoting Gauss elimination.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Uses intrinsic abs(), int() - GPU compatible.
!   - CRITICAL: Vertical data dependency in forward elimination (k-loop).
!   - Each (i,j) column is independent, but k iterations are sequential.
!   - Backward substitution also has vertical dependency.
!   - Partial pivoting version (impopt=2) has additional index indirection.
! Next:
!   - Use batched tridiagonal solver (cuSPARSE gtsv2StridedBatch).
!   - Or implement custom Thomas algorithm kernel per (i,j) column.
!   - Each column can be solved independently - batch across (i,j).
!   - Consider cyclic reduction for better parallelism if needed.
! Runtime:
!   - Calls: 14418
!   - AvgLoops: 802.8K
!   - TotalTime: 136.847s (4.59%)
!   - AvgTime: 9.491ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
