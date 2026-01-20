# Kernel 338: subroutine

## Source Location
- **File**: Src/turbuvw.f90
- **Subroutine**: subroutine
- **Line**: ~243

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates velocity turbulent mixing for u, v, w components

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 804.6K
- **Total Time**: 13.603s
- **Average Time per Call**: 37.786ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: turbuvw.f90 :: subroutine s_turbuvw
! Summary : Calculates velocity turbulent mixing for u, v, w components
!           using stress tensor divergence with terrain and map factors.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Very large parallel region with many conditional branches.
!   - Multiple code paths based on trnopt, mfcopt, mpopt, advopt.
!   - Stencil operations on stress tensors (t11, t22, t33, t12, t13, t23).
!   - Temporary arrays reused (tmp1, t11, t22 as scratch).
!   - All grid points independent within each loop nest.
! Next:
!   - Consider separating into multiple kernels by component (u, v, w).
!   - Evaluate conditions outside kernel to select specific code path.
!   - OpenACC kernels with collapse(2) on i,j loops.
!   - Data region should cover all stress tensors and force arrays.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 804.6K
!   - TotalTime: 13.603s (0.46%)
!   - AvgTime: 37.786ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
