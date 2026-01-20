# Kernel 136: s_getrich

## Source Location
- **File**: Src/getrich.f90
- **Subroutine**: s_getrich
- **Line**: ~130

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates bulk Richardson number on surface for stability

## Runtime Profile (from test_real)
- **Calls**: 377
- **Average Loop Length**: 806.4K
- **Total Time**: 0.028s
- **Average Time per Call**: 0.073ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getrich.f90 :: s_getrich
! Summary : Calculates bulk Richardson number on surface for stability
!           assessment, with special handling for sea ice regions
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic max() function - GPU compatible
!   - 2D loop over surface grid points (i,j)
!   - Conditional update for sea ice (land=1) with weighted average
!   - Module constants g, icz0m, icz0h, rchmin used from m_comphy
!   - No loop-carried dependencies
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Direct port to OpenACC parallel loop
!   - Collapse j,i loops for better GPU occupancy
!   - Ensure module constants are accessible on device
! Runtime:
!   - Calls: 377
!   - AvgLoops: 806.4K
!   - TotalTime: 0.028s (0.00%)
!   - AvgTime: 0.073ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
