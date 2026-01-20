# Kernel 180: s_instvel

## Source Location
- **File**: Src/instvel.f90
- **Subroutine**: s_instvel
- **Line**: ~142

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate maximum instantaneous wind velocity from u, v, w velocity

## Runtime Profile
- **Status**: Not executed or no runtime data available

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: instvel.f90 :: s_instvel
! Summary : Calculate maximum instantaneous wind velocity from u, v, w velocity
!           components and turbulent kinetic energy (tke)
! GPU diff: Easy
! Findings:
!   - Simple 3D loop with straightforward calculations
!   - Calls intrinsic max and sqrt functions (GPU-compatible)
!   - Private variables: k, i, j, u8s2, v8s2, w8s2
!   - Reads from u, v, w, tke arrays
!   - Updates maxvl array using max function (reduction-like pattern)
!   - No sync constructs or complex data dependencies
! Next:
!   - Can be directly ported to GPU with OpenACC parallel loop
!   - The max operation is thread-safe for independent grid points
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
