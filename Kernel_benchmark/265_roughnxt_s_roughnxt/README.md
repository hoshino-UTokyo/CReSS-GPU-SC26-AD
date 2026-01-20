# Kernel 265: s_roughnxt

## Source Location
- **File**: Src/roughnxt.f90
- **Subroutine**: s_roughnxt
- **Line**: ~116

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Update sea surface roughness length for next time step based on friction velocity

## Runtime Profile (from test_real)
- **Calls**: 361
- **Average Loop Length**: 806.4K
- **Total Time**: 0.012s
- **Average Time per Call**: 0.033ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: roughnxt.f90 :: s_roughnxt
! Summary : Update sea surface roughness length for next time step based on friction velocity
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses intrinsic max function
!   - Simple 2D loop with conditional on land use (land < 3)
!   - Conditional on ust threshold for different roughness formulas
!   - Writes to z0m and z0h arrays
!   - No synchronization constructs
! Next:
!   - Convert to OpenACC with teams distribute parallel for
!   - Straightforward GPU port with collapse(2) clause
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.012s (0.00%)
!   - AvgTime: 0.033ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
