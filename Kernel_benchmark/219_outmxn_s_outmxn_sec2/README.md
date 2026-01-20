# Kernel 219: s_outmxn

## Source Location
- **File**: Src/outmxn.f90
- **Subroutine**: s_outmxn
- **Line**: ~291
- **Section**: 2 of 2 in this subroutine

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Find grid indices (i,j,k) of max/min values using reductions

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: outmxn.f90 :: s_outmxn
! Summary : Find grid indices (i,j,k) of max/min values using reductions
!           after max/min values have been computed in previous region.
! GPU diff: Medium
! Findings:
!   - Uses reduction(max/min) for index arrays maxi,maxj,maxk,mini,minj,mink
!   - Comparison uses tolerance chkeps for floating point matching
!   - Depends on maxeps, mineps computed in previous parallel region
!   - Triple nested loop over full 3D domain
! Next:
!   - Can be fused with value-finding kernel using atomic argmax/argmin
!   - Or use two-pass: first find values, then find indices
!   - Consider storing linear index then decomposing to i,j,k
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
