# Kernel 143: s_getxy

## Source Location
- **File**: Src/getxy.f90
- **Subroutine**: s_getxy
- **Line**: ~166

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate x and y coordinates for different grid staggering

## Runtime Profile (from test_real)
- **Calls**: 2
- **Average Loop Length**: 901
- **Total Time**: 0.000s
- **Average Time per Call**: 0.113ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getxy.f90 :: s_getxy
! Summary : Calculate x and y coordinates for different grid staggering
!           (scalar, u, v, w points) based on grid spacing and MPI domain.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic real function (GPU compatible)
!   - Multiple conditional branches for different stagger types
!   - Separate 1D loops for x and y arrays
!   - No global writes, only output arrays x and y are modified
! Next:
!   - Direct translation to OpenACC with teams distribute
!   - 1D arrays are small, consider keeping on CPU or async transfer
!   - Separate kernels for x and y may be more efficient
! Runtime:
!   - Calls: 2
!   - AvgLoops: 901
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.113ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
