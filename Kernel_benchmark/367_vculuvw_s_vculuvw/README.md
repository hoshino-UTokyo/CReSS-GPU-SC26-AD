# Kernel 367: s_vculuvw

## Source Location
- **File**: Src/vculuvw.f90
- **Subroutine**: s_vculuvw
- **Line**: ~192

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Compute vertical velocity (u,v,w) advection using Cubic Lagrange scheme

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: vculuvw.f90 :: s_vculuvw
! Summary : Compute vertical velocity (u,v,w) advection using Cubic Lagrange scheme
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Reads from up,vp,wp,wc,wc8s; writes to uf,vf,wf (all 3D arrays)
!   - Three separate sections for u, v, w advection
!   - Outer k loops are serial with nested !$omp do for i,j
!   - wc8u/wc8v computed as local averages (private temporaries)
!   - Conditional branches for upwind/downwind stencil selection
! Next:
!   - Collapse k,j,i loops for each velocity component
!   - Use OpenACC teams distribute parallel do collapse(3)
!   - Map all velocity arrays to device with proper in/out semantics
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
