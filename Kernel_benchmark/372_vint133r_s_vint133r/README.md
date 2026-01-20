# Kernel 372: s_vint133r

## Source Location
- **File**: Src/vint133r.f90
- **Subroutine**: s_vint133r
- **Line**: ~144

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Interpolate variable to model grid with undefined value handling (radar data)

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vint133r.f90 :: s_vint133r
! Summary : Interpolate variable to model grid with undefined value handling (radar data)
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - getindx() called before parallel region (safe)
!   - Uses lim35n, lim34n from m_commath for undefined value markers
!   - Reads from zph, invar (3D), z1d (1D); writes to outvar (3D)
!   - First section: fill undefined values outside flat plane range
!   - Second section: interpolate with validity check on input values
!   - Conditional branches for level selection and validity checking
! Next:
!   - Collapse k,j,i loops for GPU parallelism
!   - Use OpenACC teams distribute parallel do collapse(3)
!   - Ensure lim35n, lim34n constants are accessible on device
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
