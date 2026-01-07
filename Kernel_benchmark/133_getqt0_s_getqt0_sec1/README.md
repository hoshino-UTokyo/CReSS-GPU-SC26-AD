# Kernel 133: s_getqt0

## Source Location
- **File**: Src/getqt0.f90
- **Subroutine**: s_getqt0
- **Line**: ~282
- **Section**: 1 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Initializes bubble-shaped tracer distribution using cosine function

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getqt0.f90 :: s_getqt0
! Summary : Initializes bubble-shaped tracer distribution using cosine function
!           based on distance from center points (qt0opt=1 or 2)
! GPU diff: Medium
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic cos(), sqrt(), real() functions - GPU compatible
!   - Multiple nested loops: iqt (bubbles), k, j, i
!   - Conditional write to qt array based on distance check (str < 1.0)
!   - Module variables xs, ys coordinates accessed
!   - ctr array populated inside parallel region
! Next:
!   - Port inner k,j,i loops to GPU, keep iqt loop on host or unroll
!   - Ensure xs, ys, ctr arrays are mapped to device
!   - Conditional writes may cause thread divergence
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
