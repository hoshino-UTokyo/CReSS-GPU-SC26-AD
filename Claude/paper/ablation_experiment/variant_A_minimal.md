# Phase 3: CPU Benchmark Extraction (Variant A: Minimal)

## Objective

Extract each OpenMP parallel section into a standalone benchmark program that:
1. Reads input data dumped from the actual simulation
2. Executes the kernel multiple times
3. Validates output against reference data
4. Reports timing statistics

## Target Kernels

The 15 kernels listed in `target_kernels.txt`. For each kernel, create a benchmark under `Kernel_benchmark/<id>_<name>/`.

## Dump Module

A shared dump module is available at `Kernel_benchmark/dump_kernel_data.f90`. It provides:
- `dump_init(kernel_name)` : Initialize dump directory
- `dump_params(filename, ...)` : Dump scalar parameters
- `dump_array_2d/3d/4d(filename, arr, ..)` : Dump arrays (binary)
- `dump_finalize()` : Print completion message

## Build & Run

- Build: `csh compile_radlib.csh solver compile.conf`
- Clean build (if preprocessor flags change): `csh compile_radlib.csh clean`
- Run simulation: `cd test_real && qsub MIYABI-solver.sh`
- Check job status: `qstat <jobid>`

## Simulation Config

The simulation uses 36 time steps (1/10 scale experiment).

---

## Benchmark Program Template

```fortran
program kernel_benchmark_<name>
  use omp_lib
  implicit none

  ! Variable declarations
  ! Benchmark control variables
  ! Timing variables

  ! 1. Read configuration
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! 2. Read parameters
  call read_parameters(trim(data_dir)//'/params.txt', ...)

  ! 3. Allocate arrays
  allocate(...)

  ! 4. Read input data
  call read_array_3d(trim(data_dir)//'/input.bin', ...)

  ! 5. Read reference output
  call read_array_3d(trim(data_dir)//'/output_ref.bin', ...)

  ! 6. Warmup iterations
  do iter = 1, warmup_iterations
    call kernel_<name>(...)
  end do

  ! 7. Benchmark iterations (timed)
  t_total = 0.0d0
  do iter = 1, num_iterations
    t_start = omp_get_wtime()
    call kernel_<name>(...)
    t_end = omp_get_wtime()
    t_total = t_total + (t_end - t_start)
  end do
  t_avg = t_total / dble(num_iterations)

  ! 8. Validate output
  call validate_output(output, output_ref, tolerance, passed)

  ! 9. Report results
  call print_results(t_avg, t_total, max_error, passed)

  ! 10. Cleanup
  deallocate(...)

contains

  subroutine kernel_<name>(...)
    ! Original OpenMP parallel section
  end subroutine

  ! I/O utilities: read_config, read_parameters, read_array_*, validate_output, print_results

end program
```

### Configuration File Format

`benchmark.conf`:
```
./data          # Data directory path
10              # Number of benchmark iterations
2               # Warmup iterations
1.0e-5          # Validation tolerance (relative error)
```

### Makefile Template

```makefile
FC = gfortran
FFLAGS = -O3 -fopenmp -march=native

BENCHMARK = kernel_benchmark

all: $(BENCHMARK)

$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

run: $(BENCHMARK)
	./$(BENCHMARK)

clean:
	rm -f $(BENCHMARK) *.o *.mod

.PHONY: all run clean
```

### Directory Structure (per kernel)

```
Kernel_benchmark/<id>_<name>/
├── kernel_benchmark.f90       # Benchmark program
├── Makefile
├── benchmark.conf
└── data/
    ├── params.txt
    ├── <input>.bin
    └── <output>_ref.bin
```

## Quality Checklist

- [ ] Benchmark compiles without errors
- [ ] Benchmark runs without crashes
- [ ] Validation passes (error count = 0)
- [ ] Timing is reasonable
- [ ] All input arrays read from files
- [ ] All output arrays validated
- [ ] Makefile works correctly

### Expected Output

```
==================================================
 Kernel Benchmark: <name>
==================================================
 Grid size: ni=   902, nj=   902, nk=   128
 OpenMP threads:      8
==================================================
 Average time:    16.234 ms
 Total time:     162.340 ms
 Min time:        15.891 ms
 Max time:        17.123 ms
--------------------------------------------------
 Max relative error:   1.2345E-07
 Tolerance:            1.0000E-05
 Error count:                   0
 Validation: PASSED
==================================================
```

### Special Cases

| Case | Handling |
|------|----------|
| Multiple OpenMP sections in one subroutine | Separate directories (`_sec1`, `_sec2`, etc.) |
| Kernels with reductions | Dump initial value, validate final result |
| Kernels calling external subroutines | Copy subroutine into `contains` section |

---

## Progress Report

When instructed to write a progress report, create a file at:
```
Claude/instructions/progress/progress_YYYY-MM-DD.md
```

Include the following sections:
- **Summary**: Brief description of what was accomplished
- **Completed Tasks**: List of completed items
- **In Progress**: Items currently being worked on
- **Blockers / Issues**: Problems encountered
- **Next Steps**: Priority list for next session (be specific: include file names, kernel IDs)
- **Statistics**: Kernels completed out of 15
