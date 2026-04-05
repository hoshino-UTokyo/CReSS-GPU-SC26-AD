# Phase 3: CPU Benchmark Extraction

[Back to Main](../main.md)

## Objective

Extract each OpenMP parallel section into a standalone benchmark program that:
1. Reads input data dumped from the actual simulation
2. Executes the kernel multiple times
3. Validates output against reference data
4. Reports timing statistics

---

## Two-Phase Process

| Phase | Description |
|-------|-------------|
| **Data Dump** | Modify source to dump input/output during simulation |
| **Benchmark Creation** | Create standalone program for each kernel |

---

## Phase 3.1: Data Dump Implementation

### Dump Module

Location: `Kernel_benchmark/dump_kernel_data.f90`

```fortran
module m_dump_kernel
  ! - dump_init(kernel_name)           : Initialize dump directory
  ! - dump_params(filename, ...)       : Dump scalar parameters
  ! - dump_array_2d(filename, arr, ..) : Dump 2D array (binary)
  ! - dump_array_3d(filename, arr, ..) : Dump 3D array (binary)
  ! - dump_array_4d(filename, arr, ..) : Dump 4D array (binary)
  ! - dump_finalize()                  : Print completion message
end module
```

### Instrumentation Pattern

#### Add Dump Control Variables

```fortran
use m_dump_kernel

integer, save :: dump_call_count = 0
integer, parameter :: DUMP_TARGET_CALL = <N>  ! From omp_profile.txt
logical, save :: dump_done = .false.
```

#### Add Dump Code Before OpenMP Section

```fortran
dump_call_count = dump_call_count + 1

if (dump_call_count == DUMP_TARGET_CALL .and. .not. dump_done) then
  call dump_init('<kernel_name>')
  call dump_params('params.txt', ni, nj, nk, ...)
  call dump_array_3d('input1.bin', array1, ...)
end if

!$omp parallel ...
```

#### Add Dump Code After OpenMP Section

```fortran
!$omp end parallel

if (dump_call_count == DUMP_TARGET_CALL .and. .not. dump_done) then
  call dump_array_3d('output_ref.bin', output_array, ...)
  call dump_finalize()
  dump_done = .true.
end if
```

### Variable Discovery Rule

**CRITICAL: Use `default(none)` to find all variables**

1. Set `!$omp parallel default(none) private(...) shared(...)`
2. Compile - errors reveal all variables used
3. Add sharing attributes iteratively
4. Analyze each variable's in/out role (see below)
5. Save complete list to `variable_list.txt`

### `variable_list.txt` Format

Tab-separated, one variable per line:

```
# name          type            rank    attributes      omp_clause      inout   condition
ni              integer         scalar  parameter       shared          in      -
nj              integer         scalar  parameter       shared          in      -
dt              real(8)         scalar  -               shared          in      -
u               real(8)         3d      allocatable     shared          in      -
uf              real(8)         3d      allocatable     shared          inout   -
wfrc            real(8)         3d      allocatable     shared          out     -
qall            real(8)         4d      allocatable     shared          in      -
qasl            real(8)         4d      allocatable     shared          inout   cphopt >= 2
nccn            real(8)         3d      pointer         shared          in      aslopt >= 1
k               integer         scalar  -               private         -       -
```

#### Column Definitions

| Column | Values | Description |
|--------|--------|-------------|
| `name` | variable name | |
| `type` | `integer`, `real(4)`, `real(8)`, `character`, `logical` | Fortran type |
| `rank` | `scalar`, `1d`, `2d`, `3d`, `4d` | Dimensionality |
| `attributes` | `parameter`, `allocatable`, `pointer`, `-` | Declaration attributes (comma-separated if multiple) |
| `omp_clause` | `shared`, `private`, `firstprivate`, `reduction(+:)`, etc. | OpenMP sharing clause |
| `inout` | `in`, `out`, `inout`, `-` | Role within the OpenMP region |
| `condition` | Fortran expression or `-` | Guard condition if variable is only accessed conditionally |

#### In/Out Determination Rules

1. **No function calls in region** (typical case):
   - `in`: appears only on the RHS of assignments or in conditionals
   - `out`: appears only on the LHS of assignments
   - `inout`: appears on both LHS and RHS
   - `-`: private/loop variables (not relevant for dump)

2. **Function/subroutine calls in region** (rare — see Phase 1 survey):
   - Check the called routine's `intent(in)`, `intent(out)`, `intent(inout)` for each argument
   - If intent is not declared, trace one level into the routine to determine actual usage
   - If the call chain is deeper than one level, mark as `inout` conservatively and add a note

3. **Private variables**: mark inout as `-` (not dumped, reconstructed locally in benchmark)

#### Condition Column Rules

The `condition` column prevents dump failures on unallocated/unassociated arrays:

- If a variable is accessed only inside `if (cphopt >= 2) then ... end if`, record `cphopt >= 2`
- If a variable is accessed only inside `if (allocated(var)) then ...`, record `allocated(var)`
- If a variable is always accessed unconditionally, record `-`
- The dump code must wrap the corresponding `dump_array_*` call in the same guard:

```fortran
if (cphopt >= 2) then
  call dump_array_3d('qasl.bin', qasl, ...)
end if
```

### What to Dump

| Category | How to Handle |
|----------|---------------|
| `in` / `inout` arrays | Dump as binary (`.bin`) before OpenMP region |
| `out` / `inout` arrays | Dump as binary (`.bin`) after OpenMP region (reference output) |
| `in` scalars | Dump to `params.txt` |
| `parameter` / module constants | Define as parameters in benchmark (no dump needed) |
| `private` variables | Not dumped (reconstructed locally in benchmark) |
| Conditional variables | Dump with guard condition from `variable_list.txt` |

---

## Phase 3.2: Benchmark Program Structure

### Program Template

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

---

## Directory Structure

```
Kernel_benchmark/
├── dump_kernel_data.f90          # Dump module (shared)
├── README.md                      # Index of all kernels
│
├── 001_inidef_s_inidef/
│   ├── README.md                  # Kernel info
│   ├── kernel_benchmark.f90       # Benchmark program
│   ├── Makefile
│   ├── benchmark.conf
│   └── data/
│       ├── params.txt
│       ├── input1.bin
│       └── output_ref.bin
│
└── ... (387 kernel directories)
```

---

## Makefile Template

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

---

## Workflow Summary

The process has 3 stages. Stages 1 and 2 are done once for all kernels; stage 3 is per-kernel.

### Stage 1: Variable List Creation (per kernel)

For each kernel, create `variable_list.txt`:

1. **Identify the OpenMP region** in the source file
2. **Use `default(none)`** to discover all variables (see Variable Discovery Rule above)
3. **Classify each variable** (type, rank, omp_clause, inout, condition)
4. **Save** as `Kernel_benchmark/<id>_<name>/variable_list.txt`

### Stage 2: Dump Data Collection (one simulation run)

Add dump instrumentation for **all target kernels at once**, then run the simulation **once** to collect all dump data:

1. **Add dump code** to each source file based on its `variable_list.txt`
2. **Build** the instrumented simulation
3. **Run simulation once** — all kernels dump their data in a single run
4. **Copy dump files** to each kernel's `data/` directory
5. **Remove dump instrumentation** from source files

### Stage 3: Benchmark Creation (per kernel)

For each kernel:

1. **Create benchmark program** following template, using `variable_list.txt` for declarations and I/O
2. **Build and test**: `make && ./kernel_benchmark`
3. **Verify validation passes** (error count = 0)

---

## Special Cases

| Case | Handling |
|------|----------|
| Multiple OpenMP sections in one subroutine | Separate directories (`_sec1`, `_sec2`, etc.) |
| Kernels with reductions | Dump initial value, validate final result |
| Kernels calling external subroutines | Copy subroutine into `contains` section |
| Unexecuted kernels (Count=0) | Use different config or create synthetic data |

---

## Quality Checklist

- [ ] Benchmark compiles without errors
- [ ] Benchmark runs without crashes
- [ ] Validation passes (error count = 0)
- [ ] Timing is reasonable
- [ ] All input arrays read from files
- [ ] All output arrays validated
- [ ] README.md updated with results
- [ ] Makefile works correctly

---

## Expected Output

```
==================================================
 Kernel Benchmark: diver3d
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

---

## Related Documents

- [Phase 2: Runtime Profiling](02_runtime_profiling.md) - Provides call counts
- [Phase 4: GPU Benchmark](04_gpu_benchmark.md) - Next step

---

*Original: benchmark_instruction.md*
