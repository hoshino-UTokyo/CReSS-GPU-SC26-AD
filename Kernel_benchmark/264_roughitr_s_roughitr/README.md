# Kernel 264: s_roughitr

## Source Location
- **File**: Src/roughitr.f90
- **Subroutine**: s_roughitr
- **Line**: ~185

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Calculate sea surface roughness length by iteration using bulk coefficient

## Runtime Profile (from test_real)
- **Calls**: 16
- **Average Loop Length**: 806.4K
- **Total Time**: 0.001s
- **Average Time per Call**: 0.033ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: roughitr.f90 :: s_roughitr
! Summary : Calculate sea surface roughness length by iteration using bulk coefficient
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - Located inside iterate: do loop (outer iteration loop is serial)
!   - Uses intrinsic functions: abs, max
!   - Conditional logic based on land use (land < 3) and ust threshold
!   - Updates z0m, z0h, dz0m arrays
!   - External subroutine calls getrich, bulksfc, chkitr before/after parallel region
!   - No synchronization constructs inside parallel region
! Next:
!   - Convert to OpenACC with teams distribute parallel for
!   - Keep iteration control on host, only offload inner 2D loop
!   - May need to manage z0m, z0h data between iterations on GPU
! Runtime:
!   - Calls: 16
!   - AvgLoops: 806.4K
!   - TotalTime: 0.001s (0.00%)
!   - AvgTime: 0.033ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
