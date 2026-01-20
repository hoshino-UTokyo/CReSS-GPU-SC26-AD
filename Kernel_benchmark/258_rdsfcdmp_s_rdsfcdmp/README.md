# Kernel 258: s_rdsfcdmp

## Source Location
- **File**: Src/rdsfcdmp.f90
- **Subroutine**: s_rdsfcdmp
- **Line**: ~570

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Set constant land use categories based on terrain height and sfcopt setting

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: rdsfcdmp.f90 :: s_rdsfcdmp
! Summary : Set constant land use categories based on terrain height and sfcopt setting
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Two separate do loops based on sfcopt condition (conditional inside parallel region)
!   - Simple element-wise assignment to land array
!   - No synchronization constructs or reductions
! Next:
!   - Convert to OpenACC with teams distribute parallel for
!   - Move sfcopt conditional outside kernel for simpler GPU code
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
