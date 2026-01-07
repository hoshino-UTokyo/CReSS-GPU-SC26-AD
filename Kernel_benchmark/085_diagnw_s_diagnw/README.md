# Kernel 085: s_diagnw

## Source Location
- **File**: Src/diagnw.f90
- **Subroutine**: s_diagnw
- **Line**: ~135

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute diagnostic concentrations of cloud water and rain water

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.4M
- **Total Time**: 2.073s
- **Average Time per Call**: 5.759ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diagnw.f90 :: s_diagnw
! Summary : Compute diagnostic concentrations of cloud water and rain water
!           based on base state density and water hydrometeor mixing ratios.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region (only intrinsics: sqrt, min, max)
!   - No global/module variable writes
!   - No synchronization constructs (barriers, critical, atomic)
!   - Simple 3D loop with k-loop outside, j-i loops inside with schedule(runtime)
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - Consider collapse(3) after loop restructuring for better GPU utilization
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 2.073s (0.07%)
!   - AvgTime: 5.759ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
