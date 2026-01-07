# Kernel 220: s_outpbl

## Source Location
- **File**: Src/outpbl.f90
- **Subroutine**: s_outpbl
- **Line**: ~332

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute 10m wind, 1.5m pressure/temperature/humidity, surface temp,

## Runtime Profile (from test_real)
- **Calls**: 4
- **Average Loop Length**: 806.4K
- **Total Time**: 0.001s
- **Average Time per Call**: 0.215ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: outpbl.f90 :: s_outpbl
! Summary : Compute 10m wind, 1.5m pressure/temperature/humidity, surface temp,
!           cloud cover, and surface fluxes for output diagnostics.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region (pure arithmetic)
!   - No global/module variable writes (only intent(inout) arrays)
!   - No sync constructs (barrier, critical, atomic)
!   - Conditional branches based on fmois (dry vs moist) with separate do loops
! Next:
!   - Convert to OpenACC or OpenACC with data directives for arrays
!   - Collapse nested i,j loops for better GPU occupancy
! Runtime:
!   - Calls: 4
!   - AvgLoops: 806.4K
!   - TotalTime: 0.001s (0.00%)
!   - AvgTime: 0.215ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
