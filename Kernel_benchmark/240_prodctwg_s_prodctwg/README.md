# Kernel 240: s_prodctwg

## Source Location
- **File**: Src/prodctwg.f90
- **Subroutine**: s_prodctwg
- **Line**: ~221

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate graupel production rate for wet/dry growth processes

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 0.986s
- **Average Time per Call**: 0.022ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: prodctwg.f90 :: s_prodctwg
! Summary : Calculate graupel production rate for wet/dry growth processes
!           in cloud microphysics scheme.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Complex conditional logic with many if/else branches
!   - Reads/writes to multiple 3D arrays (clir, clis, clig, clsr, clsg, etc.)
!   - Many private local variables (pgdry, lfice, cligw, clsgw, sink, a)
!   - Data-dependent branching may cause thread divergence on GPU
!   - No explicit barriers but implicit at !$omp end do
! Next:
!   - Use OpenACC kernels or parallel loop with private clause for locals
!   - Thread divergence from conditionals may reduce GPU efficiency
!   - Consider restructuring conditionals to minimize divergence
!   - Ensure all read/write arrays are in data region
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 0.986s (0.03%)
!   - AvgTime: 0.022ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
