# Kernel 289: s_setname

## Source Location
- **File**: Src/setname.f90
- **Subroutine**: s_setname
- **Line**: ~546

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Copies land use parameters from arrays to name table indexed by category,

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 10
- **Total Time**: 0.000s
- **Average Time per Call**: 0.122ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: setname.f90 :: s_setname
! Summary : Copies land use parameters from arrays to name table indexed by category,
!           including albedo, beta, roughness, heat capacity, and diffusivity.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple 1D loops copying from land category arrays to name table
!   - All loops are independent with private iid index
!   - No synchronization constructs
!   - Loop count is numctg_lnd (small, number of land categories)
! Next:
!   - Straightforward GPU port with simple 1D kernel
!   - Consider keeping on CPU due to small loop count
!   - If porting, use single kernel for all table copies
! Runtime:
!   - Calls: 1
!   - AvgLoops: 10
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.122ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
