# Kernel 198: s_lsps0

## Source Location
- **File**: Src/lsps0.f90
- **Subroutine**: s_lsps0
- **Line**: ~151

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Apply lateral sponge damping to optional scalar variable forcing term

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: lsps0.f90 :: s_lsps0
! Summary : Apply lateral sponge damping to optional scalar variable forcing term
!           with optional smoothing, damping to initial state
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to sfrc (output forcing term) and tmp1 (temporary array)
!   - Multiple worksharing constructs with branching logic for lspopt
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Collapse nested i,j loops for better GPU occupancy
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
