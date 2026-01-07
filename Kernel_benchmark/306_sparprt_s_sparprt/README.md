# Kernel 306: s_sparprt

## Source Location
- **File**: Src/sparprt.f90
- **Subroutine**: s_sparprt
- **Line**: ~149

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Separates pressure and potential temperature into base state and

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: sparprt.f90 :: s_sparprt
! Summary : Separates pressure and potential temperature into base state and
!           perturbation values by subtracting interpolated base state
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Single !$omp do loop with schedule(runtime)
!   - Simple element-wise subtraction operations
!   - Writes to ppdat, ptpdat arrays (in-place update)
!   - No synchronization constructs besides implicit barriers
! Next:
!   - Data managed automatically via Unified Memory
!   - Convert to !$acc parallel loop collapse(3)
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
