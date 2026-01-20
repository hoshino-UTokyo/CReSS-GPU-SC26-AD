# Kernel 066: s_copy2d

## Source Location
- **File**: Src/copy2d.f90
- **Subroutine**: s_copy2d
- **Line**: ~108

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Simple 2D array copy from invar to outvar.

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: copy2d.f90 :: s_copy2d
! Summary : Simple 2D array copy from invar to outvar.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls
!   - Trivial memory copy operation
!   - Independent element-wise operations
! Next:
!   - Straightforward GPU port with loop collapse
!   - Consider using device-to-device memcpy for efficiency
!   - May be better to keep data resident on GPU and avoid copy calls
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
