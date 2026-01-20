# Kernel 203: s_mk2cgbw

## Source Location
- **File**: Src/mk2cgbw.f90
- **Subroutine**: s_mk2cgbw
- **Line**: ~130

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Convert units from [m][k] to [c][g] for warm bin cloud physics arrays

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: mk2cgbw.f90 :: s_mk2cgbw
! Summary : Convert units from [m][k] to [c][g] for warm bin cloud physics arrays
!           (density, pressure, mixing ratios, concentrations, precipitation)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple scalar multiplication operations
!   - Writes to rbr, rst, rbv, p, qwbin, nwbin, prr arrays (in-place modification)
!   - Nested loops over bin categories (n) and vertical levels (k)
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Collapse nested i,j loops for better GPU occupancy
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
