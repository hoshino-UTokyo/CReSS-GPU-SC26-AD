# Kernel 233: s_phvbcuvw

## Source Location
- **File**: Src/phvbcuvw.f90
- **Subroutine**: s_phvbcuvw
- **Line**: ~397

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Calculate differential phase speed terms for u,v,w velocity

## Runtime Profile (from test_real)
- **Calls**: 360
- **Average Loop Length**: 1
- **Total Time**: 10.499s
- **Average Time per Call**: 29.163ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: phvbcuvw.f90 :: s_phvbcuvw
! Summary : Calculate differential phase speed terms for u,v,w velocity
!           components at open boundary conditions (west/east/south/north).
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
!   - Consider using OpenACC kernels with private scalars bc0, bc1, bc2
!   - Boundary-only computation may benefit from separate small kernels
! Runtime:
!   - Calls: 360
!   - AvgLoops: 1
!   - TotalTime: 10.499s (0.35%)
!   - AvgTime: 29.163ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
