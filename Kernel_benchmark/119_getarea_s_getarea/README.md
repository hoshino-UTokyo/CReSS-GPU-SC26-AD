# Kernel 119: s_getarea

## Source Location
- **File**: Src/getarea.f90
- **Subroutine**: s_getarea
- **Line**: ~234

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Compute total area of each boundary plane (top/bottom and

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 806.4K
- **Total Time**: 0.000s
- **Average Time per Call**: 0.120ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getarea.f90 :: s_getarea
! Summary : Compute total area of each boundary plane (top/bottom and
!           lateral) by summing grid cell areas with map scale factors.
! GPU diff: Medium
! Findings:
!   - Multiple reduction operations (+ for area0, areaw, areae, areas, arean)
!   - Complex conditional logic based on mfcopt, mpopt, boundary flags
!   - Different loops for different boundary planes
!   - Uses MPI module variables (ebw, ebe, ebs, ebn, isub, jsub, etc.)
! Next:
!   - Use GPU reduction kernels for area summation
!   - May need separate kernels for each boundary plane
!   - Consider whether GPU overhead is worthwhile for boundary-only computation
! Runtime:
!   - Calls: 1
!   - AvgLoops: 806.4K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.120ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
