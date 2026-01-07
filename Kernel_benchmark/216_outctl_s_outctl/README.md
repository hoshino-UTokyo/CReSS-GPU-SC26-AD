# Kernel 216: s_outctl

## Source Location
- **File**: Src/outctl.f90
- **Subroutine**: s_outctl
- **Line**: ~350

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate constant height z1d at scalar points for GrADS

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: outctl.f90 :: s_outctl
! Summary : Calculate constant height z1d at scalar points for GrADS
!           control file output, interpolating stretched coordinates.
! GPU diff: Easy
! Findings:
!   - Small loop over k from 2 to nk-2
!   - Conditionals on fproc, mype, dmplev checked outside omp do
!   - Simple arithmetic: z1d(k) = 0.5*(zsth(k)+zsth(k+1))
!   - Only executed on root process (mype.eq.root)
!   - No inter-thread dependencies
! Next:
!   - Very small computation, likely not worth GPU offload
!   - Keep on CPU as it runs only on root process
!   - Simple loop could be left as serial if needed
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
