# Kernel 179: s_initund

## Source Location
- **File**: Src/initund.f90
- **Subroutine**: s_initund
- **Line**: ~235

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Initialize soil and sea temperature arrays (tund, tundp) based on

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 806.4K
- **Total Time**: 0.003s
- **Average Time per Call**: 3.217ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: initund.f90 :: s_initund
! Summary : Initialize soil and sea temperature arrays (tund, tundp) based on
!           surface type (land/sea), SST data, and atmospheric conditions
! GPU diff: Medium
! Findings:
!   - Multiple conditional branches (sfcopt, sfcdat, advopt, land type)
!   - Calls intrinsic exp, log, min, real functions (GPU-compatible)
!   - Private variables: k, i, j
!   - Reads from land, sst, pbr, ptbr, pp, ptp, ek arrays
!   - Writes to tund and tundp arrays
!   - Uses shared ek array computed within parallel region
!   - Multiple omp do regions within single parallel block
!   - No sync constructs between threads
! Next:
!   - Can be ported to GPU with OpenACC parallel loop
!   - Ensure ek array is properly handled (computed then used)
!   - Consider separating different sfcopt cases into different kernels
! Runtime:
!   - Calls: 1
!   - AvgLoops: 806.4K
!   - TotalTime: 0.003s (0.00%)
!   - AvgTime: 3.217ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
