# Kernel 231: s_pgradiv

## Source Location
- **File**: Src/pgradiv.f90
- **Subroutine**: s_pgradiv
- **Line**: ~147

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Compute vertical pressure gradient force using implicit method,

## Runtime Profile (from test_real)
- **Calls**: 14400
- **Average Loop Length**: 100.4M
- **Total Time**: 117.847s
- **Average Time per Call**: 8.184ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: pgradiv.f90 :: s_pgradiv
! Summary : Compute vertical pressure gradient force using implicit method,
!           updating forcing term for w equation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region (pure arithmetic)
!   - No global/module variable writes (only intent(inout) fw, fpdvj arrays)
!   - No sync constructs
!   - Two separate k-loops: first computes fpdvj, second updates fw
!   - Second loop has k-1 dependency on fpdvj (read from previous level)
! Next:
!   - Convert to OpenACC with collapse for k,j,i loops
!   - First loop is independent; second loop needs fpdvj from k-1 level
!   - Can fuse loops or ensure proper synchronization between them
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 100.4M
!   - TotalTime: 117.847s (3.95%)
!   - AvgTime: 8.184ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
