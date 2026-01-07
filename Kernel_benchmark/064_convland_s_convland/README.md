# Kernel 064: s_convland

## Source Location
- **File**: Src/convland.f90
- **Subroutine**: s_convland
- **Line**: ~104

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Converts land use categories between integer and real types,

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: convland.f90 :: s_convland
! Summary : Converts land use categories between integer and real types,
!           either real(land)+0.1 or nint(rland).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls (only intrinsic nint, real)
!   - Simple 2D loop with straightforward type conversion
!   - Two mutually exclusive branches based on fproc string
!   - Independent grid point operations
! Next:
!   - Straightforward GPU port with OpenACC or OpenACC
!   - Collapse i,j loops for better occupancy
!   - Consider data movement optimization if called frequently
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
