# Kernel 023: s_allocuni

## Source Location
- **File**: Src/allocuni.f90
- **Subroutine**: s_allocuni
- **Line**: ~226

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Initialize arrays for the unite program (tmp1-4, iodmp, var)

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: allocuni.f90 :: s_allocuni
! Summary : Initialize arrays for the unite program (tmp1-4, iodmp, var)
!           to zero for file merging operations
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to module-level arrays tmp1-4, iodmp, var from m_comuni
!   - Multiple simple initialization loops with no data dependencies
!   - Four separate do loops for different array dimensions
! Next:
!   - Straightforward GPU port with OpenACC parallel loops
!   - Collapse nested loops in var initialization
!   - Small arrays (nk, nio_uni) may not benefit from GPU offload
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
