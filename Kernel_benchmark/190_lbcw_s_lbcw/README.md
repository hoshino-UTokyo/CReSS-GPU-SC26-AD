# Kernel 190: s_lbcw

## Source Location
- **File**: Src/lbcw.f90
- **Subroutine**: s_lbcw
- **Line**: ~161

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Apply lateral boundary conditions for z-component velocity (w)

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 113.1K
- **Total Time**: 0.058s
- **Average Time per Call**: 0.004ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: lbcw.f90 :: s_lbcw
! Summary : Apply lateral boundary conditions for z-component velocity (w)
!           by copying values from interior to boundary points.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls within parallel region
!   - Simple array copy operations (wf(boundary) = wf(interior))
!   - Multiple conditionals based on boundary type (wbc,ebc,sbc,nbc) and advopt
!   - Uses module variables from m_commpi (ebw,ebe,ebs,ebn,isub,jsub,nisub,njsub)
!   - No synchronization constructs besides implicit barrier at omp end do
! Next:
!   - Convert to OpenACC or OpenACC kernels
!   - Consider collapsing k and j/i loops for better GPU utilization
!   - Data managed automatically via Unified Memory
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 113.1K
!   - TotalTime: 0.058s (0.00%)
!   - AvgTime: 0.004ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
