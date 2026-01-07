# Kernel 384: s_vsps0

## Source Location
- **File**: Src/vsps0.f90
- **Subroutine**: s_vsps0
- **Line**: ~118

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Applies vertical sponge damping to optional scalar forcing

## Runtime Profile (from test_real)
- **Calls**: 1440
- **Average Loop Length**: 101.2M
- **Total Time**: 7.870s
- **Average Time per Call**: 5.466ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vsps0.f90 :: s_vsps0
! Summary : Applies vertical sponge damping to optional scalar forcing
!           term, relaxing variable toward zero (initial state).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No conditional branches inside parallel region
!   - k loop starts from ksp0(2)-1 (variable start index)
!   - Simple arithmetic update to sfrc array
!   - No synchronization constructs other than implicit barriers
! Next:
!   - Straightforward GPU port with collapse on j,i loops
!   - Handle variable k-range start with appropriate kernel bounds
!   - Map sfrc, sp, rbct, rst arrays to device
! Runtime:
!   - Calls: 1440
!   - AvgLoops: 101.2M
!   - TotalTime: 7.870s (0.26%)
!   - AvgTime: 5.466ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
