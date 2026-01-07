# Kernel 315: subroutine

## Source Location
- **File**: Src/stepwi.f90
- **Subroutine**: subroutine
- **Line**: ~323

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Prepares coefficient matrices for vertical implicit solver

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 99.5M
- **Total Time**: 417.326s
- **Average Time per Call**: 28.981ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: stepwi.f90 :: subroutine s_stepwi
! Summary : Prepares coefficient matrices for vertical implicit solver
!           of w-equation, computing tridiagonal matrix elements.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region (pure arithmetic only).
!   - Reads module constant g from comphy.
!   - No synchronization constructs.
!   - Conditional on buyopt selects different physics formulation.
!   - Sets up tridiagonal system (tmp1=lower, tmp2=diag, tmp3=upper).
!   - Note: gaussel/gseidel solvers called outside parallel region.
! Next:
!   - Coefficient setup is embarrassingly parallel - easy to port.
!   - Tridiagonal solver (gaussel) needs separate GPU implementation
!     (batched tridiagonal solver or cyclic reduction).
!   - Consider cuSPARSE gtsv2 or custom kernel for vertical solve.
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 99.5M
!   - TotalTime: 417.326s (14.01%)
!   - AvgTime: 28.981ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
