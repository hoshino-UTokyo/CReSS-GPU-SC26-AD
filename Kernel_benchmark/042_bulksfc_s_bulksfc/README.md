# Kernel 042: s_bulksfc

## Source Location
- **File**: Src/bulksfc.f90
- **Subroutine**: s_bulksfc
- **Line**: ~217

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Calculates bulk coefficients for surface momentum and heat

## Runtime Profile (from test_real)
- **Calls**: 385
- **Average Loop Length**: 806.4K
- **Total Time**: 0.289s
- **Average Time per Call**: 0.750ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: bulksfc.f90 :: s_bulksfc
! Summary : Calculates bulk coefficients for surface momentum and heat
!           fluxes using Monin-Obukhov similarity theory with stability.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Heavy use of intrinsics: sqrt, log, exp, acos, atan, cos, tan, abs, min, max
!   - Complex nested conditionals (land type, stability, ice coverage)
!   - Significant control flow divergence based on rch sign and land values
!   - Many local temporary variables (a through f, cmice, chice, dz0m, dz0h)
!   - Single 2D loop over surface grid points
! Next:
!   - Convert to OpenACC with collapse(2) for i,j loops
!   - Map all 2D input/output arrays to device
!   - GPU divergence may reduce efficiency; consider separating cases
!   - Transcendental functions may benefit from fast-math approximations
! Runtime:
!   - Calls: 385
!   - AvgLoops: 806.4K
!   - TotalTime: 0.289s (0.01%)
!   - AvgTime: 0.750ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
