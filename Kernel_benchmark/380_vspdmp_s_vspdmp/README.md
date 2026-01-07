# Kernel 380: s_vspdmp

## Source Location
- **File**: Src/vspdmp.f90
- **Subroutine**: s_vspdmp
- **Line**: ~180

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Calculates relaxed vertical sponge damping coefficients with

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 126
- **Total Time**: 0.010s
- **Average Time per Call**: 10.285ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vspdmp.f90 :: s_vspdmp
! Summary : Calculates relaxed vertical sponge damping coefficients with
!           maximum z-coordinate search and cosine-based damping profiles.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Uses max() intrinsic for reduction-like operation on z1dmax
!   - !$omp single block for sequential ksp0 index search (do_k_1, do_k_2)
!   - Multiple k loops with different purposes (init, max-find, coef-calc)
!   - Conditional vspopt branches inside parallel region
!   - Potential race condition in z1dmax(k)=max(...) without proper reduction
! Next:
!   - z1dmax computation needs reduction or atomic operations for GPU
!   - Sequential ksp0 search should remain on CPU or use parallel reduction
!   - Split into separate kernels: max-find, ksp0-search, coef-calculation
!   - Use cosine from device math library
! Runtime:
!   - Calls: 1
!   - AvgLoops: 126
!   - TotalTime: 0.010s (0.00%)
!   - AvgTime: 10.285ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
