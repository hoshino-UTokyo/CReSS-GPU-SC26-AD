# Kernel 112: s_forcept

## Source Location
- **File**: Src/forcept.f90
- **Subroutine**: s_forcept
- **Line**: ~326

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute potential temperature pt by adding base state ptbr

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.4M
- **Total Time**: 1.361s
- **Average Time per Call**: 3.780ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: forcept.f90 :: s_forcept
! Summary : Compute potential temperature pt by adding base state ptbr
!           and perturbation ptpp for turbulent mixing calculation.
! GPU diff: Easy
! Findings:
!   - Simple element-wise addition of two arrays
!   - No function calls inside parallel region
!   - No reductions or synchronization
!   - Only writes to pt array
! Next:
!   - Direct GPU kernel port with straightforward 3D mapping
!   - Consider fusing with subsequent turbulent mixing kernels
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 1.361s (0.05%)
!   - AvgTime: 3.780ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
