# Kernel 294: s_sheartke

## Source Location
- **File**: Src/sheartke.f90
- **Subroutine**: s_sheartke
- **Line**: ~116

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates shear production term in TKE equation by adding

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 100.4M
- **Total Time**: 1.960s
- **Average Time per Call**: 5.446ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: sheartke.f90 :: s_sheartke
! Summary : Calculates shear production term in TKE equation by adding
!           Jacobian times eddy viscosity times deformation squared to forcing.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple 3D loop with single multiply-add operation per point
!   - All loops independent with private i,j,k indices
!   - No synchronization constructs
!   - Minimal computation per grid point
! Next:
!   - Straightforward GPU port with 3D kernel
!   - Use OpenACC/OpenACC with collapse(3)
!   - Good candidate for kernel fusion with other TKE terms
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.4M
!   - TotalTime: 1.960s (0.07%)
!   - AvgTime: 5.446ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
