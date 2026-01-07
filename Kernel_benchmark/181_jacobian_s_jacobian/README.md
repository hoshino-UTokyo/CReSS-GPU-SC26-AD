# Kernel 181: s_jacobian

## Source Location
- **File**: Src/jacobian.f90
- **Subroutine**: s_jacobian
- **Line**: ~223

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate transformation Jacobian components (j31, j32, jcb) from

## Runtime Profile (from test_real)
- **Calls**: 1
- **Average Loop Length**: 103.6M
- **Total Time**: 0.008s
- **Average Time per Call**: 8.288ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: jacobian.f90 :: s_jacobian
! Summary : Calculate transformation Jacobian components (j31, j32, jcb) from
!           physical coordinates (x, y, z, zph)
! GPU diff: Easy
! Findings:
!   - Three separate loop nests computing j31, j32, and jcb arrays
!   - Simple arithmetic operations (subtraction, division)
!   - Private variables: k, i, j
!   - Reads from x, y, z, zph arrays
!   - Writes to j31, j32, jcb arrays
!   - No function calls within the parallel region
!   - No sync constructs or data dependencies between grid points
! Next:
!   - Can be directly ported to GPU with OpenACC parallel loop
!   - Consider fusing loops for better GPU memory access patterns
! Runtime:
!   - Calls: 1
!   - AvgLoops: 103.6M
!   - TotalTime: 0.008s (0.00%)
!   - AvgTime: 8.288ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
