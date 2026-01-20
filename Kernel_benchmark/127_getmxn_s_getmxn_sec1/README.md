# Kernel 127: s_getmxn

## Source Location
- **File**: Src/getmxn.f90
- **Subroutine**: s_getmxn
- **Line**: ~200
- **Section**: 1 of 4 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Finds maximum value across MPI gathered buffer and identifies

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: getmxn.f90 :: s_getmxn
! Summary : Finds maximum value across MPI gathered buffer and identifies
!           the corresponding index location
! GPU diff: Medium
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic abs(), max(), sign() functions - GPU compatible
!   - Contains reduction operations (max: maxvl, maxeps)
!   - Iterates over MPI processor elements (npe), typically small count
!   - Accesses module buffers mxnbuf, idxbuf from m_combuf
!   - Module constant eps used from m_commath
! Next:
!   - Small loop size (npe) may not benefit from GPU offloading
!   - Consider keeping on CPU or using atomic operations
!   - Reduction clause supported in OpenACC
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
