# Kernel 082: s_diagnci

## Source Location
- **File**: Src/diagnci.f90
- **Subroutine**: s_diagnci
- **Line**: ~126

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate diagnostic cloud ice concentrations from ice mixing

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diagnci.f90 :: s_diagnci
! Summary : Calculate diagnostic cloud ice concentrations from ice mixing
!           ratio using a simple linear scaling with inverse max mass.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Private variable k for outer loop
!   - Writes to nidia output array
!   - Simple point-wise multiplication operation
!   - Completely independent iterations
! Next:
!   - Direct conversion to OpenACC with collapsed loops
!   - Minimal data transfer: input qice, output nidia
!   - Excellent GPU candidate due to simplicity
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
