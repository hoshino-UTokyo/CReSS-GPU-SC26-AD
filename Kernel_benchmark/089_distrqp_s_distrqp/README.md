# Kernel 089: s_distrqp

## Source Location
- **File**: Src/distrqp.f90
- **Subroutine**: s_distrqp
- **Line**: ~176

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Distribute radar-observed precipitation mixing ratios to rain,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: distrqp.f90 :: s_distrqp
! Summary : Distribute radar-observed precipitation mixing ratios to rain,
!           snow, graupel, and hail categories based on model hydrometeor ratios.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region (only intrinsic: abs used before region)
!   - No global/module variable writes
!   - No synchronization constructs
!   - Multiple branches (datype_rdr, cphopt, haiopt) but all are data-parallel
!   - Conditional distribution based on thresholds
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - Conditionals per grid point are fine for GPU
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
