# Kernel 073: s_curveuvw

## Source Location
- **File**: Src/curveuvw.f90
- **Subroutine**: s_curveuvw
- **Line**: ~175

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate earth curvature forcing terms for u, v, w velocity

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 102.4M
- **Total Time**: 17.817s
- **Average Time per Call**: 49.492ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: curveuvw.f90 :: s_curveuvw
! Summary : Calculate earth curvature forcing terms for u, v, w velocity
!           equations using temporary arrays and map scale factors.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Private variable k used for outer loop
!   - Writes to ufrc, vfrc, wfrc (forcing terms) and tmp1-tmp5 (temporaries)
!   - Multiple sequential k-loops with data dependencies between them
!   - Conditional branches based on mpopt and mfcopt options
! Next:
!   - Convert to OpenACC with data region for all arrays
!   - May need to fuse some k-loops or restructure for better parallelism
!   - Handle conditional logic for map projection options on GPU
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 17.817s (0.60%)
!   - AvgTime: 49.492ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
