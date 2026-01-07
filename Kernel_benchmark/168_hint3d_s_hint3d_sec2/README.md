# Kernel 168: s_hint3d

## Source Location
- **File**: Src/hint3d.f90
- **Subroutine**: s_hint3d
- **Line**: ~386
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Perform 3D horizontal interpolation (linear or parabolic) from

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: hint3d.f90 :: s_hint3d
! Summary : Perform 3D horizontal interpolation (linear or parabolic) from
!           data grid to model grid for each vertical level.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region; uses intrinsics only
!   - Branching based on intopt (linear vs parabolic) and mpopt
!   - Parabolic interpolation has more complex stencil access
!   - Outer k-loop with inner j,i loops parallelized via !$omp do
!   - Writes to var 3D array (output)
!   - No sync constructs
! Next:
!   - Convert to OpenACC with collapse(3) on k,j,i loops
!   - Consider separate kernels for linear vs parabolic interpolation
!   - Parabolic case may benefit from shared memory for stencil
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
