# Kernel 224: s_pbldrv

## Source Location
- **File**: Src/pbldrv.f90
- **Subroutine**: s_pbldrv
- **Line**: ~243

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Convert virtual potential temperature back to potential temperature

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: pbldrv.f90 :: s_pbldrv
! Summary : Convert virtual potential temperature back to potential temperature
!           perturbation after PBL diffusion calculations.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region (pure arithmetic)
!   - No global/module variable writes (only intent(inout) ptp array)
!   - No sync constructs
!   - Conditional branches based on fmois (dry vs moist) with different formulas
!   - Loop over k levels with nested i,j loops
! Next:
!   - Convert to OpenACC with collapse for k,j,i loops
!   - Can be combined with preceding PBL subroutine calls into single kernel
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
