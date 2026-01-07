# Kernel 102: diver3d (s_diver3d)

## Source Location
- **File**: Src/diver3d.f90
- **Subroutine**: s_diver3d
- **Line**: ~184

## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Calculate 3D negative divergence using velocity components weighted by Jacobian and map scale factors.

## Runtime Profile (from test_real)
- **Calls**: 14,400
- **Average Loop Length**: 102.5M iterations
- **Total Time**: 232.922s (7.82% of total)
- **Average Time per Call**: 16.175ms

## Quick Start

```bash
# Build
make

# Generate test data (902x902x128 grid, ~4GB)
make data

# Run benchmark
make run

# Or run with specific thread count
OMP_NUM_THREADS=8 ./kernel_benchmark
```

## Files

| File | Description |
|------|-------------|
| `kernel_benchmark.f90` | Main benchmark program |
| `generate_data.f90` | Test data generator |
| `Makefile` | Build configuration |
| `benchmark.conf` | Benchmark configuration (auto-generated) |
| `data/` | Input/output data directory (auto-generated) |

## Benchmark Structure

```
kernel_benchmark.f90
├── program kernel_benchmark_diver3d
│   ├── read_config()        - Load benchmark settings
│   ├── read_parameters()    - Load grid/physics params
│   ├── read_array_*()       - Load input arrays
│   ├── warmup loop          - CPU cache warmup
│   ├── benchmark loop       - Timed kernel execution
│   ├── validation           - Compare with reference
│   └── report               - Print results
└── contains
    └── kernel_diver3d()     - The actual kernel (OpenMP)
```

## Configuration

Edit `benchmark.conf` to customize:
```
./data          # Data directory
10              # Number of benchmark iterations
2               # Warmup iterations
1.0e-5          # Validation tolerance
```

Or modify `generate_data.f90` to change grid size:
```fortran
ni = 902    ! X dimension
nj = 902    ! Y dimension
nk = 128    ! Z dimension (vertical levels)
```

## Output Example

```
==================================================
 Kernel Benchmark: diver3d
==================================================
 Grid size: ni=  902, nj=  902, nk=  128
 Options: mpopt=    10, mfcopt=     1
 Warmup iterations:      2
 Benchmark iterations:     10
 OpenMP threads:      8
==================================================
 Loading input data...
 Loading reference output...
 Running warmup iterations...
 Running benchmark iterations...
 Validating output...

==================================================
 Results
==================================================
 Average time:    15.234567 ms
 Total time:     152.345670 ms
 Min time:        14.891234 ms
 Max time:        16.123456 ms
--------------------------------------------------
 Max relative error:   1.2345E-07
 Tolerance:            1.0000E-05
 Error count:                   0
 Validation: PASSED
==================================================
```

## Thread Scaling Test

```bash
make run-scaling
```

This runs the benchmark with 1, 2, 4, 8, 16, 32 threads.

## Notes

- Grid size 902×902×128 matches the actual simulation configuration
- Data files are ~4GB total (binary format for efficiency)
- Reference output is computed serially for correctness
- Validation uses relative error with configurable tolerance
