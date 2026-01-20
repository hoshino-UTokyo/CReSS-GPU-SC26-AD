# Kernel 047: s_cg2mkbw

## Source Location
- **File**: Src/cg2mkbw.f90
- **Subroutine**: s_cg2mkbw
- **Line**: ~121

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Converts units from CGS [c][g] to MKS [m][k] for warm bin

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: cg2mkbw.f90 :: s_cg2mkbw
! Summary : Converts units from CGS [c][g] to MKS [m][k] for warm bin
!           cloud microphysics arrays (mass, concentration, precipitation).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple scalar multiplication operations
!   - Nested loops over bin categories (n), vertical levels (k), and grid (i,j)
!   - Independent updates to mwbin, nwbin, prr arrays
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Use collapse for nested loops (consider collapse(3) or collapse(4))
!   - May need to handle n-loop separately if nqw varies at runtime
!   - Simple structure well-suited for GPU offload
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
