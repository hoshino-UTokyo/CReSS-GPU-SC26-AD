# Kernel 093: s_diverpih

## Source Location
- **File**: Src/diverpih.f90
- **Subroutine**: s_diverpih
- **Line**: ~202

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate horizontal divergence for pressure equation (HEVI method),

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 100.4M
- **Total Time**: 48.758s
- **Average Time per Call**: 3.386ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: diverpih.f90 :: s_diverpih
! Summary : Calculate horizontal divergence for pressure equation (HEVI method),
!           including terrain-following coordinate contributions when trnopt>0.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - External call to diver2d before this parallel region (already annotated)
!   - No global/module variable writes
!   - No synchronization constructs
!   - Multiple branches (trnopt, mfcopt, mpopt) with different loop structures
!   - Flat terrain case is simple; curved grid has multi-stage computation
!   - Uses pdiv(i,j,nk) as temporary storage in one branch
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - tmp1, tmp2, tmp3 temporaries need GPU allocation
!   - diver2d call should also be GPU-ported for full offload
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 100.4M
!   - TotalTime: 48.758s (1.64%)
!   - AvgTime: 3.386ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
