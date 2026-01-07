# Kernel 309: s_steppts

## Source Location
- **File**: Src/steppts.f90
- **Subroutine**: s_steppts
- **Line**: ~216

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Advances potential temperature perturbation in time using forcing

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: steppts.f90 :: s_steppts
! Summary : Advances potential temperature perturbation in time using forcing
!           and gravity wave terms
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Single !$omp do loop with schedule(runtime)
!   - Simple element-wise update: ptpf += dts*(ptfrc+ptsml)/rst
!   - Writes only to ptpf array
!   - No synchronization constructs besides implicit barriers
! Next:
!   - Data managed automatically via Unified Memory
!   - Convert to !$acc parallel loop collapse(2)
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
