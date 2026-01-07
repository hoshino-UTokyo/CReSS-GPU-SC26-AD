# Kernel 078: s_depsit

## Source Location
- **File**: Src/depsit.f90
- **Subroutine**: s_depsit
- **Line**: ~245

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate evaporation rate from rain to vapor and deposition

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 7.399s
- **Average Time per Call**: 0.162ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: depsit.f90 :: s_depsit
! Summary : Calculate evaporation rate from rain to vapor and deposition
!           rates from vapor to ice hydrometeors (cloud ice, snow, graupel).
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic exp, int, log, max functions (GPU-compatible)
!   - Accesses lookup tables ckoe, pkoe indexed by temperature
!   - Private variable k for outer loop; many private local scalars
!   - Writes to vdvr, vdvi, vdvs, vdvg output arrays
!   - Complex conditional logic based on temperature and mixing ratios
!   - Special case handling for nk=1 vs nk>1
! Next:
!   - Convert to OpenACC with data region for all input/output arrays
!   - Copy lookup tables ckoe, pkoe to device
!   - May need to restructure conditionals for GPU efficiency
!   - Consider separating nk=1 case into distinct kernel
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 7.399s (0.25%)
!   - AvgTime: 0.162ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
