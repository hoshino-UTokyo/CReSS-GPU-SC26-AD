# Kernel 382: s_vspqv

## Source Location
- **File**: Src/vspqv.f90
- **Subroutine**: s_vspqv
- **Line**: ~163

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Applies vertical sponge damping to water vapor mixing ratio

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 101.2M
- **Total Time**: 2.742s
- **Average Time per Call**: 7.616ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vspqv.f90 :: s_vspqv
! Summary : Applies vertical sponge damping to water vapor mixing ratio
!           forcing term, relaxing to GPV data or base state value.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Conditional branch (vspopt, gpvvar) selects damping target
!   - k loop starts from ksp0-1 (variable start index)
!   - Simple arithmetic update to qvfrc array
!   - No synchronization constructs other than implicit barriers
! Next:
!   - Straightforward GPU port with collapse on j,i loops
!   - Handle variable k-range start with appropriate kernel bounds
!   - Map qvfrc, qvp, qvgpv, qvtd, qvbr, rbct, rst arrays to device
! Runtime:
!   - Calls: 360
!   - AvgLoops: 101.2M
!   - TotalTime: 2.742s (0.09%)
!   - AvgTime: 7.616ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
