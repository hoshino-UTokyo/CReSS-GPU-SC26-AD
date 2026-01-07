# Kernel 222: s_paractl

## Source Location
- **File**: Src/paractl.f90
- **Subroutine**: s_paractl
- **Line**: ~483
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Find minimum/maximum latitude and longitude from corner points

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: paractl.f90 :: s_paractl
! Summary : Find minimum/maximum latitude and longitude from corner points
!           using reduction operations.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region (pure arithmetic with min/max)
!   - No global/module variable writes
!   - Uses reduction(min:) and reduction(max:) for latmin, lonmin, latmax, lonmax
!   - Small loop iteration count (4-5 iterations) - may not benefit from GPU
! Next:
!   - Consider keeping on CPU due to small iteration count
!   - If porting, use GPU reduction primitives or atomic operations
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
