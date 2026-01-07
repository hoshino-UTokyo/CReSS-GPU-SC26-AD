# Kernel 204: s_more0q

## Source Location
- **File**: Src/more0q.f90
- **Subroutine**: s_more0q
- **Line**: ~243

## Analysis
- **GPU Difficulty**: Medium
- **Summary**: Force mixing ratios to be non-negative by scaling microphysical process

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 806.4K
- **Total Time**: 2.650s
- **Average Time per Call**: 0.058ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: more0q.f90 :: s_more0q
! Summary : Force mixing ratios to be non-negative by scaling microphysical process
!           rates when sink exceeds available mass for cloud water, ice, rain, snow, graupel
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic max() function - GPU compatible
!   - Uses module variable evapor_opt from m_temparam
!   - Complex conditional logic with multiple branching per hydrometeor type
!   - Writes to many microphysical rate arrays (nuvi, nuci, clcr, clcs, ... etc)
!   - Branching on nk==1 for 2D vs 3D handling
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - May need to restructure conditionals for better GPU branch divergence
!   - Collapse nested i,j loops for better GPU occupancy
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 2.650s (0.09%)
!   - AvgTime: 0.058ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
