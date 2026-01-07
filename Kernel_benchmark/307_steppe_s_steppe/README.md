# Kernel 307: s_steppe

## Source Location
- **File**: Src/steppe.f90
- **Subroutine**: s_steppe
- **Line**: ~222

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Advances pressure perturbation in time using forcing and acoustic

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: steppe.f90 :: s_steppe
! Summary : Advances pressure perturbation in time using forcing and acoustic
!           terms with horizontally explicit/vertically explicit method
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Single !$omp do loop with schedule(runtime)
!   - Simple element-wise update: ppf += dts*(pfrc+psml)/jcb
!   - Writes only to ppf array
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
