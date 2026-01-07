# Kernel 057: s_cnt2phy

## Source Location
- **File**: Src/cnt2phy.f90
- **Subroutine**: s_cnt2phy
- **Line**: ~166

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Convert contravariant vertical velocity (wc) to physical vertical

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: cnt2phy.f90 :: s_cnt2phy
! Summary : Convert contravariant vertical velocity (wc) to physical vertical
!           velocity (w) accounting for terrain and map scale factors
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function/subroutine calls inside parallel region
!   - No reductions or synchronization constructs
!   - Multiple conditional branches based on trnopt, sthopt, mfcopt, mpopt
!   - Intermediate arrays j31u2, j32v2, mf25 computed and used within region
!   - Simple element-wise arithmetic operations
! Next:
!   - Direct OpenACC with collapse(2) for inner loops
!   - Consider fusing loops where possible to reduce kernel launches
!   - Data dependencies between j31u2/j32v2 computation and w computation
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
