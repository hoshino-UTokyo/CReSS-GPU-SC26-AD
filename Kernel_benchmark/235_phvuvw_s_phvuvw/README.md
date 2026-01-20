# Kernel 235: s_phvuvw

## Source Location
- **File**: Src/phvuvw.f90
- **Subroutine**: s_phvuvw
- **Line**: ~361

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Calculate velocity phase speed (u,v,w) for open boundary

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 112.2K
- **Total Time**: 5.198s
- **Average Time per Call**: 14.440ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: phvuvw.f90 :: s_phvuvw
! Summary : Calculate velocity phase speed (u,v,w) for open boundary
!           conditions on all four boundaries (west, east, south, north).
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region (only intrinsic abs, sign, min, max, mod)
!   - Writes to output arrays ucpx, ucpy, vcpx, vcpy, wcpx, wcpy
!   - Uses shared work arrays cpavex, cpavey, u8v, v8u for vertical averaging
!   - Complex conditional branching based on boundary condition options
!   - Multiple sequential k-loops with workshared inner j/i loops
!   - No explicit barriers but implicit at !$omp end do
! Next:
!   - Collapse nested loops where possible for better GPU occupancy
!   - Consider using OpenACC kernels directive with appropriate private clauses
!   - Boundary-only computation may benefit from separate small kernels
! Runtime:
!   - Calls: 360
!   - AvgLoops: 112.2K
!   - TotalTime: 5.198s (0.17%)
!   - AvgTime: 14.440ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
