# Kernel 280: s_setblk

## Source Location
- **File**: Src/setblk.f90
- **Subroutine**: s_setblk
- **Line**: ~309

## Analysis
- **GPU Difficulty**: Hard
- **Summary**: Calculate thermodynamic properties (T, saturation, latent heat),

## Runtime Profile (from test_real)
- **Calls**: 45720
- **Average Loop Length**: 1
- **Total Time**: 26.486s
- **Average Time per Call**: 0.579ms

## Original Annotation
```fortran
!@llm start meta_info ----------------------------------------------------
! Location: setblk.f90 :: s_setblk
! Summary : Calculate thermodynamic properties (T, saturation, latent heat),
!           air properties (viscosity, conductivity), and hydrometeor parameters
! GPU diff: Hard
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic functions (exp, log, sqrt) extensively
!   - Complex conditional structure: nk==1 vs nk>1, abs(cphopt)<=3 vs ==4
!   - Many output arrays (t, tcel, qvsw, qvsi, lv, ls, lf, kp, mu, dv, etc.)
!   - Uses 2D work array nu(0:ni+1,0:nj+1) between loop nests within k loop
!   - Uses module constants from m_commath, m_comphy
!   - Multiple !$omp do regions within single parallel region
! Next:
!   - Consider restructuring to avoid nu dependency between loop nests
!   - Large number of output arrays requires careful data management
!   - Branch structure may benefit from separate GPU kernels per case
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 1
!   - TotalTime: 26.486s (0.89%)
!   - AvgTime: 0.579ms
!@llm end meta_info ------------------------------------------------------
```

## Benchmark Files (TODO)
- `kernel.f90` - Extracted kernel code
- `driver.f90` - Benchmark driver
- `data/` - Input data for benchmark
- `Makefile` - Build configuration
