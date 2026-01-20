# Kernel 310: s_steps

## Source Location
- **File**: Src/steps.f90
- **Subroutine**: s_steps
- **Line**: ~518

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Advances all scalar variables (ptp, qv, hydrometeors, aerosols,

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 100.4M
- **Total Time**: 16.586s
- **Average Time per Call**: 46.072ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: steps.f90 :: s_steps
! Summary : Advances all scalar variables (ptp, qv, hydrometeors, aerosols,
!           tracers, TKE) to next time step using forcing terms
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Many !$omp do loops with schedule(runtime)
!   - Complex conditional branching (fmois, cphopt, haiopt, qcgopt, aslopt, etc.)
!   - Writes to dtdrst, ptpf, qvf, qwtrf, nwtrf, qicef, nicef, qcwtrf, qcicef, qaslf, qtf, tkef
!   - Uses max/min intrinsics for clipping values
!   - n_sub loop variable for array dimension iteration
!   - No synchronization constructs besides implicit barriers
! Next:
!   - Data managed automatically via Unified Memory
!   - Consider separating each variable update into distinct kernels
!   - Branching may require conditional kernel launches or unified kernels
!   - Use collapse(2) for nested loops
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.4M
!   - TotalTime: 16.586s (0.56%)
!   - AvgTime: 46.072ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
