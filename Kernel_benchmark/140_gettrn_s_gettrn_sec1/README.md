# Kernel 140: s_gettrn

## Source Location
- **File**: Src/gettrn.f90
- **Subroutine**: s_gettrn
- **Line**: ~189
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Initialize terrain height array with constant flat value

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 810.0K
- **Total Time**: 0.000s
- **Average Time per Call**: 0.015ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: gettrn.f90 :: s_gettrn (trnopt=0 branch)
! Summary : Initialize terrain height array with constant flat value
!           using max of mountain height+base or sea surface height.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic max function (GPU compatible)
!   - Simple constant assignment to all grid points
!   - No global writes, only output array ht is modified
! Next:
!   - Direct translation to OpenACC with teams distribute
!   - Consider using GPU memset-like operation for constant fill
! Runtime:
!   - Calls: 1
!   - AvgLoops: 810.0K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.015ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
