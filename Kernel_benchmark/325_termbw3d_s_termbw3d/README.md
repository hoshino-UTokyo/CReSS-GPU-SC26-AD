# Kernel 325: s_termbw3d

## Source Location
- **File**: Src/termbw3d.f90
- **Subroutine**: s_termbw3d
- **Line**: ~184

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate terminal velocity for water bins using empirical

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: termbw3d.f90 :: s_termbw3d
! Summary : Calculate terminal velocity for water bins using empirical
!           formulas based on droplet size regimes (3D version)
! GPU diff: Medium
! Findings:
!   - Three size regimes with different formulas (r<1e-3, r<5.35e-2, r<0.35)
!   - Multiple exp/log intrinsic calls per grid point
!   - Coefficients c1-c5 computed only when ncp==1 (first bin)
!   - Conditional branches based on droplet radius
!   - 3D arrays with additional k dimension vs 2D version
! Next:
!   - Compute c1-c5 coefficients once on first call (ncp==1)
!   - Consider separating size regimes to reduce divergence
!   - Collapse k,j,i loops for GPU parallelization
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
