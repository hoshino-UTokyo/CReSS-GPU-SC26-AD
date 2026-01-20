# Kernel 295: s_shedding

## Source Location
- **File**: Src/shedding.f90
- **Subroutine**: s_shedding
- **Line**: ~142

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculates shedding rates of liquid water from snow and graupel to rain,

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 1.185s
- **Average Time per Call**: 0.026ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: shedding.f90 :: s_shedding
! Summary : Calculates shedding rates of liquid water from snow and graupel to rain,
!           based on temperature and collection/production rates.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses module constant t0 from m_comphy for temperature threshold
!   - Conditional branches based on mixing ratio thresholds and temperature
!   - Handles nk=1 case separately (2D) vs nk>1 case (3D)
!   - All loops independent with private i,j,k indices
!   - No synchronization constructs
! Next:
!   - Straightforward GPU port with conditional logic preserved
!   - Use OpenACC/OpenACC with collapse for nested loops
!   - Consider single kernel handling both nk cases with runtime check
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.185s (0.04%)
!   - AvgTime: 0.026ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
