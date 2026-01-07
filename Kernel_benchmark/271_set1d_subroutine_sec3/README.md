# Kernel 271: subroutine

## Source Location
- **File**: Src/set1d.f90
- **Subroutine**: subroutine
- **Line**: ~449
- **Section**: 3 of 8 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculates virtual temperature and integrates Exner function

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: set1d.f90 :: subroutine s_set1d (Exner/virtual temp calc)
! Summary : Calculates virtual temperature and integrates Exner function
!           hydrostatically for pressure-based sounding data.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constant epsav from comphy.
!   - Uses !$omp single for sequential integration (vertical dependency).
!   - Uses reduction(min:) for pimin error checking.
! Next:
!   - Keep sequential integration on host (small nlev).
!   - Similar approach as temperature-based sounding.
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
