# Kernel 182: s_kh8uv

## Source Location
- **File**: Src/kh8uv.f90
- **Subroutine**: s_kh8uv
- **Line**: ~138

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Set horizontal eddy diffusivity at u and v points (rkh8u, rkh8v)

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.3M
- **Total Time**: 3.068s
- **Average Time per Call**: 8.522ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: kh8uv.f90 :: s_kh8uv
! Summary : Set horizontal eddy diffusivity at u and v points (rkh8u, rkh8v)
!           by averaging rkh values, with map scale factor corrections
! GPU diff: Medium
! Findings:
!   - Multiple conditional branches (mfcopt, mpopt)
!   - Simple arithmetic operations (addition, multiplication)
!   - Private variables: k, i, j
!   - Reads from rkh, rmf arrays
!   - Writes to rkh, rkh8u, rkh8v arrays
!   - rkh is modified in-place then used (potential ordering concern)
!   - Multiple omp do regions within single parallel block
!   - No sync constructs between threads
! Next:
!   - Can be ported to GPU with OpenACC parallel loop
!   - Careful attention needed for rkh modification ordering
!   - Consider separating different mpopt/mfcopt cases into different kernels
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.3M
!   - TotalTime: 3.068s (0.10%)
!   - AvgTime: 8.522ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
