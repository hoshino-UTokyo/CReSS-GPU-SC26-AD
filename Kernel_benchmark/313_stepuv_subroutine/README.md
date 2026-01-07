# Kernel 313: subroutine

## Source Location
- **File**: Src/stepuv.f90
- **Subroutine**: subroutine
- **Line**: ~292

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Time integration of u and v velocity components using

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 100.5M
- **Total Time**: 155.820s
- **Average Time per Call**: 10.821ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: stepuv.f90 :: subroutine s_stepuv
! Summary : Time integration of u and v velocity components using
!           forcing terms and acoustic mode contributions.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Simple arithmetic: uf = uf + dts * (ufrc + usml) / rst8u.
!   - All grid points are independent (embarrassingly parallel).
!   - Note: Multiple subroutine calls before/after for boundary conditions
!     and MPI communication - those need separate GPU handling.
! Next:
!   - Direct OpenACC kernels for the time stepping loops.
!   - Consider fusing u and v updates into single kernel.
!   - MPI communication calls outside parallel region need GPU-aware MPI.
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 100.5M
!   - TotalTime: 155.820s (5.23%)
!   - AvgTime: 10.821ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
