# Kernel 217: s_outdmp

## Source Location
- **File**: Src/outdmp.f90
- **Subroutine**: s_outdmp
- **Line**: ~342

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate constant height z1d array for dump output based

## Runtime Profile (from test_real)
- **Calls**: 4
- **Average Loop Length**: 125
- **Total Time**: 0.000s
- **Average Time per Call**: 0.019ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: outdmp.f90 :: s_outdmp
! Summary : Calculate constant height z1d array for dump output based
!           on dmplev option (uniform dz spacing or stretched coords).
! GPU diff: Easy
! Findings:
!   - Small loop over k from 2 to nk-2
!   - Conditionals on fdmp and dmplev checked outside omp do
!   - Simple arithmetic: z1d(k) = dz*(real(k)-1.5) or zsth average
!   - No inter-thread dependencies; each k independent
!   - Typically small nk dimension
! Next:
!   - May not benefit from GPU due to small loop size
!   - Direct port is straightforward if needed
!   - Consider keeping on CPU for simplicity
! Runtime:
!   - Calls: 4
!   - AvgLoops: 125
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.019ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
