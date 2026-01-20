# Kernel 041: s_bruntv

## Source Location
- **File**: Src/bruntv.f90
- **Subroutine**: s_bruntv
- **Line**: ~205

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates Brunt-Vaisala frequency squared for atmospheric

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.4M
- **Total Time**: 12.063s
- **Average Time per Call**: 33.508ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bruntv.f90 :: s_bruntv
! Summary : Calculates Brunt-Vaisala frequency squared for atmospheric
!           stability using potential temperature and moisture fields.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Intrinsic functions used: exp, log (transcendental, may need GPU libs)
!   - Multiple conditional branches based on fmois and cphopt (divergence)
!   - Multiple work arrays (pt, ptv, a, t) modified in sequence
!   - Some data dependencies between loop nests (e.g., pt used to compute ptv)
!   - 2D temporary array t(i,j) reused across k iterations
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Use collapse(2) for nested i,j loops
!   - May need to restructure k-loop to avoid thread-local t array issues
!   - Consider separating dry/moist cases into different GPU kernels
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 12.063s (0.40%)
!   - AvgTime: 33.508ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
