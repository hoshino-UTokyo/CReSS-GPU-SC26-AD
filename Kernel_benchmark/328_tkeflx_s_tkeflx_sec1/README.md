# Kernel 328: s_tkeflx

## Source Location
- **File**: Src/tkeflx.f90
- **Subroutine**: s_tkeflx
- **Line**: ~211
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate j31*tke and j32*tke products for terrain-following

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: tkeflx.f90 :: s_tkeflx (terrain preprocessing)
! Summary : Calculate j31*tke and j32*tke products for terrain-following
!           coordinate transformation
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls
!   - Writes to j31tke and j32tke arrays
!   - Simple stencil computations
!   - No synchronization constructs within parallel region
! Next:
!   - Straightforward GPU port
!   - Can be fused with main flux calculation if data dependencies allow
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
