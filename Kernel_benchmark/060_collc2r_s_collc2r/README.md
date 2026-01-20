# Kernel 060: s_collc2r

## Source Location
- **File**: Src/collc2r.f90
- **Subroutine**: s_collc2r
- **Line**: ~144

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates collection rate between cloud water and rain water

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: collc2r.f90 :: s_collc2r
! Summary : Calculates collection rate between cloud water and rain water
!           using exponential/logarithmic formulas for microphysics.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region (only intrinsic exp, log)
!   - Simple 3D loop with private loop variables and local temporaries
!   - No synchronization constructs beyond implicit barrier at end do
!   - Read/write to qcf and qrf arrays with independent grid points
! Next:
!   - Can be ported directly with OpenACC or OpenACC parallel loop
!   - Consider collapsing the k,j,i loops for better GPU occupancy
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
