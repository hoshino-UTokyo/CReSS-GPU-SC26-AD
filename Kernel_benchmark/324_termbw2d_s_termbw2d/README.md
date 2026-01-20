# Kernel 324: s_termbw2d

## Source Location
- **File**: Src/termbw2d.f90
- **Subroutine**: s_termbw2d
- **Line**: ~177

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate terminal velocity for water bins using empirical

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: termbw2d.f90 :: s_termbw2d
! Summary : Calculate terminal velocity for water bins using empirical
!           formulas based on droplet size regimes (2D version)
! GPU diff: Medium
! Findings:
!   - Three size regimes with different formulas (r<1e-3, r<5.35e-2, r<0.35)
!   - Multiple exp/log intrinsic calls per grid point
!   - Precomputed coefficients c1-c5 reused across bins
!   - Conditional branches based on droplet radius
! Next:
!   - Compute c1-c5 coefficients once, then loop over bins
!   - Consider separating size regimes to reduce divergence
!   - Keep coefficient arrays on GPU for reuse
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
