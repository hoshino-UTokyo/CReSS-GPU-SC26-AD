# Kernel 374: s_vint31g

## Source Location
- **File**: Src/vint31g.f90
- **Subroutine**: s_vint31g
- **Line**: ~172

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Vertically interpolates 3D GPV data to 1D flat plane with

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vint31g.f90 :: s_vint31g
! Summary : Vertically interpolates 3D GPV data to 1D flat plane with
!           extrapolation handling based on surface reference options.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Multiple conditional branches (refsfc_gpv, etrvar_gpv) control flow
!   - Nested k/kd loops with !$omp do on inner jd,id loops
!   - No synchronization constructs other than implicit barriers
!   - Output array varef written conditionally based on z-level comparisons
! Next:
!   - Collapse nested loops where possible for better GPU occupancy
!   - Consider restructuring conditionals outside parallel region
!   - Use data directives for varef, zdat, vardat, zlow arrays
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
