# Kernel 074: s_dbz2kg

## Source Location
- **File**: Src/dbz2kg.f90
- **Subroutine**: s_dbz2kg
- **Line**: ~123

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Convert radar reflectivity from dBZe to precipitation mixing

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: dbz2kg.f90 :: s_dbz2kg
! Summary : Convert radar reflectivity from dBZe to precipitation mixing
!           ratio in kg/m^3 using exponential/logarithmic transformation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic exp and log functions (GPU-compatible)
!   - Simple 3D loop with independent point-wise operations
!   - Writes to qpdat array (in-place modification)
!   - Conditional check on lim34n threshold
! Next:
!   - Direct conversion to OpenACC or OpenACC with collapsed loops
!   - Data managed automatically via Unified Memory
!   - No synchronization needed between iterations
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
