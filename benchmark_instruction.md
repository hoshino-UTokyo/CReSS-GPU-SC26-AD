# Kernel Benchmark Extraction Guide (Claude Code Instructions)

## Purpose

Extract each OpenMP parallel section into a standalone benchmark that:
1. **Reads input data** dumped from the actual simulation
2. **Executes the kernel** a specified number of times
3. **Validates output** against reference data from the simulation
4. **Reports timing statistics** for performance analysis

This enables isolated testing and GPU porting of individual kernels.

---

## Overview

### Total Kernels to Extract
- **387 OpenMP parallel sections** across 338 source files
- Profiling data available in `test_real/omp_profile.txt`
- Kernel directories prepared in `Kernel_benchmark/`

### Two-Phase Process

**Phase 1: Data Dump**
- Modify source files to dump input/output data during simulation
- Run simulation once to generate data files
- Data is dumped at the **final call** of each kernel (using call count from profiling)

**Phase 2: Benchmark Creation**
- Create standalone benchmark program for each kernel
- Single-file Fortran program with kernel code in `contains` section
- Reads dumped data, executes kernel, validates results

---

## Phase 1: Data Dump Implementation

### 1.1 Dump Module

Use the dump module located at `Kernel_benchmark/dump_kernel_data.f90`:

```fortran
module m_dump_kernel
  ! Provides:
  ! - dump_init(kernel_name)           : Initialize dump directory
  ! - dump_params(filename, ...)       : Dump scalar parameters
  ! - dump_array_2d(filename, arr, ..) : Dump 2D array (binary)
  ! - dump_array_3d(filename, arr, ..) : Dump 3D array (binary)
  ! - dump_array_4d(filename, arr, ..) : Dump 4D array (binary)
  ! - dump_finalize()                  : Print completion message
end module
```

### 1.2 Instrumentation Pattern

For each OpenMP parallel section, add dump code that triggers on the **final call**.

#### Step 1: Add module reference
```fortran
use m_dump_kernel
```

#### Step 2: Add dump control variables
```fortran
! Dump control (in variable declaration section)
integer, save :: dump_call_count = 0
integer, parameter :: DUMP_TARGET_CALL = <N>  ! From omp_profile.txt "Count" column
logical, save :: dump_done = .false.
```

#### Step 3: Add dump code before OpenMP section
```fortran
! Increment call counter
dump_call_count = dump_call_count + 1

! Dump on final call
if (dump_call_count == DUMP_TARGET_CALL .and. .not. dump_done) then
  call dump_init('<kernel_name>')
  call dump_params('params.txt', ni, nj, nk, ...)
  call dump_array_3d('input1.bin', array1, ...)
  call dump_array_3d('input2.bin', array2, ...)
  ! ... dump all input arrays
end if

!$omp parallel ...
```

#### Step 4: Add dump code after OpenMP section
```fortran
!$omp end parallel

! Dump output on final call
if (dump_call_count == DUMP_TARGET_CALL .and. .not. dump_done) then
  call dump_array_3d('output_ref.bin', output_array, ...)
  call dump_finalize()
  dump_done = .true.
end if
```

### 1.3 Variable Discovery Rule (default(none) Method)

**CRITICAL: Follow this procedure to ensure no variables are missed when dumping.**

#### Step 1: Set default(none) on the target OpenMP section
```fortran
!$omp parallel default(none) private(...) shared(...)
```

#### Step 2: Compile and collect errors
The compiler will report all variables that need explicit sharing clauses:
```
Error: Symbol 'xxx' is used but not explicitly specified in DEFAULT(NONE) clause
```

#### Step 3: Add sharing attributes iteratively
- Compile → Fix one error → Compile again → Repeat
- This systematically reveals ALL variables used in the parallel region

#### Step 4: Create variable list
Save the complete variable list to `Kernel_benchmark/<kernel_dir>/variable_list.txt`:
```
# Variables for <kernel_name>
# Discovered via default(none) method

## Shared variables (inputs/outputs)
var1  - description
var2  - description

## Private variables
i, j, k - loop indices

## Reduction variables
sum_val - reduction(+)
```

#### Step 5: Repeat for all OpenMP parallel sections
Apply this procedure to every `!$omp parallel` block in the target subroutine.

**This method guarantees complete variable coverage and prevents dump omissions.**

---

### 1.4 Determining DUMP_TARGET_CALL

Get the call count from `test_real/omp_profile.txt`:

```
   ID  File              Subroutine         Count    ...
  102  diver3d.f90       s_diver3d          14400    ...
```

For `diver3d`, set `DUMP_TARGET_CALL = 14400`.

**Important**: For kernels with Count=0 (unexecuted), skip data dump or use a different test configuration that exercises that code path.

### 1.4 What to Dump

**CRITICAL RULE: Dump ALL variables used inside the OpenMP parallel section.**

Every variable that appears inside `!$omp parallel ... !$omp end parallel` must be either:
1. Dumped to a file and read in the benchmark, OR
2. Computed from dumped parameters using the exact same formula

#### Categories of Variables to Dump

##### 1. Arrays (dump as binary files)
- All arrays read by the kernel (intent(in), intent(inout))
- All work/temporary arrays that contain meaningful input values
- Arrays computed by called subroutines (e.g., if kernel calls `diver3d` to compute `tmp1`, dump `tmp1`)

##### 2. Scalar Parameters (dump to params.txt)
- Grid dimensions (ni, nj, nk)
- Physics options (mpopt, mfcopt, trnopt, divopt, isoopt, etc.)
- Physical constants from namelist (dx, dy, dz, dxiv, dyiv, dziv, dts, etc.)

##### 3. Derived Constants (dump OR document computation)
- Constants computed from parameters: `dz05 = 0.5 * dz`, `ds308 = 0.125 * dx * dy * dz`
- Mathematical constants from modules: `oned3 = 1.0/3.0`, `eps = 1.0e-20`
- **Option A**: Dump these values to params.txt
- **Option B**: Document the formula in benchmark and recompute (must match exactly)

##### 4. Module Constants (MUST handle explicitly)
- Constants from `m_commath`: `oned3`, `eps`, `cc`, etc.
- Constants from `m_comphy`: physical constants
- **These MUST be either dumped or defined as parameters in the benchmark**

#### Output Data (dump AFTER kernel execution)
- All arrays written by the kernel (intent(out), intent(inout))
- This becomes the reference for validation

#### Common Mistakes to Avoid

1. **Missing work arrays**: If a kernel uses `tmp1` computed by another subroutine (e.g., `diver3d`), you MUST dump `tmp1` after that subroutine runs
2. **Missing module constants**: Constants like `oned3`, `eps` from `m_commath` must be defined in the benchmark
3. **Missing derived parameters**: Values like `ds308 = 0.125*dx*dy*dz` must be either dumped or computed identically
4. **Assuming zero initialization**: Work arrays may contain non-zero values from previous computations

#### Parameter File Format (params.txt)
```
<ni>           ! Grid dimension X
<nj>           ! Grid dimension Y
<nk>           ! Grid dimension Z
<option1>      ! Physics option 1
<option2>      ! Physics option 2
<real_param1>  ! Real parameter (ES20.12 format)
<real_param2>  ! Real parameter
...
```

#### Array File Format (*.bin)
- Binary stream format (`access='stream', form='unformatted'`)
- Fortran column-major order
- Single precision (real) or as declared in source

---

## Phase 2: Benchmark Program Structure

### 2.1 Single-File Program Template

Each benchmark is a self-contained Fortran program:

```fortran
!***********************************************************************
! Kernel Benchmark: <kernel_name>
!***********************************************************************
! Source: Src/<filename>.f90
! Subroutine: <subroutine_name>
! Profile ID: <id>
!***********************************************************************
program kernel_benchmark_<name>
  use omp_lib
  implicit none

  ! === Variable declarations ===
  ! (All arrays and parameters needed by the kernel)

  ! === Benchmark control ===
  integer :: num_iterations, warmup_iterations
  real :: tolerance
  character(len=256) :: data_dir

  ! === Timing ===
  real(8) :: t_start, t_end, t_total, t_avg

  ! === Main program flow ===

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

  !-----------------------------------------------------------------
  ! The kernel (extracted from original source)
  !-----------------------------------------------------------------
  subroutine kernel_<name>(...)
    ! === Original OpenMP parallel section ===
    !$omp parallel ...
    ...
    !$omp end parallel
  end subroutine

  !-----------------------------------------------------------------
  ! If kernel calls other subroutines, include them here
  !-----------------------------------------------------------------
  subroutine helper_function(...)
    ...
  end subroutine

  !-----------------------------------------------------------------
  ! I/O utilities
  !-----------------------------------------------------------------
  subroutine read_config(...)
  subroutine read_parameters(...)
  subroutine read_array_2d(...)
  subroutine read_array_3d(...)
  subroutine validate_output(...)
  subroutine print_results(...)

end program
```

### 2.2 Handling Subroutine Calls Inside Kernels

If the OpenMP parallel section calls other subroutines:

1. **Check if the subroutine is called inside `!$omp parallel`**
2. **Include the subroutine in `contains` section**
3. **Recursively include any subroutines it calls**

Example:
```fortran
contains

  subroutine kernel_pgrad(...)
    !$omp parallel ...
    call diver3d(...)  ! Called inside parallel region
    !$omp end parallel
  end subroutine

  ! Must include diver3d
  subroutine diver3d(...)
    ...
  end subroutine

end program
```

### 2.3 Benchmark Configuration File

`benchmark.conf` format:
```
./data          # Data directory path
10              # Number of benchmark iterations
2               # Warmup iterations
1.0e-5          # Validation tolerance (relative error)
```

### 2.4 Validation Approach

```fortran
subroutine validate_output(output, reference, tol, passed, max_err, err_count)
  real, intent(in) :: output(:,:,:), reference(:,:,:)
  real, intent(in) :: tol
  logical, intent(out) :: passed
  real, intent(out) :: max_err
  integer, intent(out) :: err_count

  real :: rel_err
  integer :: i, j, k

  max_err = 0.0
  err_count = 0

  do k = 1, size(output,3)
    do j = 1, size(output,2)
      do i = 1, size(output,1)
        if (abs(reference(i,j,k)) > 1.0e-30) then
          rel_err = abs(output(i,j,k) - reference(i,j,k)) / abs(reference(i,j,k))
        else
          rel_err = abs(output(i,j,k) - reference(i,j,k))
        end if
        if (rel_err > max_err) max_err = rel_err
        if (rel_err > tol) err_count = err_count + 1
      end do
    end do
  end do

  passed = (err_count == 0)
end subroutine
```

---

## Directory Structure

```
Kernel_benchmark/
├── dump_kernel_data.f90          # Dump module (shared)
├── README.md                      # Index of all kernels
│
├── 001_inidef_s_inidef/
│   ├── README.md                  # Kernel info and instructions
│   ├── kernel_benchmark.f90       # Standalone benchmark program
│   ├── Makefile                   # Build configuration
│   ├── benchmark.conf             # Runtime configuration
│   └── data/                      # Dumped data files
│       ├── params.txt
│       ├── input1.bin
│       ├── input2.bin
│       └── output_ref.bin
│
├── 102_diver3d_s_diver3d/         # Example: completed
│   ├── kernel_benchmark.f90
│   ├── generate_data.f90          # (Optional: synthetic data generator)
│   ├── Makefile
│   └── data/
│
└── ... (387 kernel directories)
```

---

## Makefile Template

```makefile
# Compiler settings
FC = gfortran
FFLAGS = -O3 -fopenmp -march=native

# Targets
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

### For Each Kernel (387 total):

1. **Identify kernel** from `Kernel_benchmark/<id>_<name>/README.md`
   - Source file location
   - Call count from profiling
   - Input/output arrays

2. **Add dump instrumentation** to source file
   - Add `use m_dump_kernel`
   - Add dump control variables with correct `DUMP_TARGET_CALL`
   - Add dump calls before/after OpenMP section

3. **Run simulation** to generate dump files
   - Execute `test_real` or appropriate test case
   - Verify dump files created in `kernel_dump/<name>/`

4. **Copy dump files** to benchmark directory
   ```bash
   cp -r kernel_dump/<name>/* Kernel_benchmark/<id>_<name>/data/
   ```

5. **Create benchmark program**
   - Extract kernel code from source
   - Create `kernel_benchmark.f90` following template
   - Include any called subroutines in `contains`

6. **Build and test**
   ```bash
   cd Kernel_benchmark/<id>_<name>/
   make
   ./kernel_benchmark
   ```

7. **Verify validation passes**
   - Max relative error should be < tolerance
   - Error count should be 0

---

## Special Cases

### Multiple OpenMP Sections in One Subroutine

Some subroutines have multiple `!$omp parallel` sections (e.g., `hculuvw` has 6).

- Each section gets its own kernel directory (`_sec1`, `_sec2`, etc.)
- Dump each section's inputs/outputs separately
- Use section-specific call counters

### Kernels with Reductions

For kernels using `!$omp reduction`:
- Dump the reduction variable initial value as input
- Dump the final reduction result as output
- Ensure benchmark initializes reduction variable before each iteration

### Kernels Calling External Subroutines

If kernel calls subroutines defined in other modules:
1. Copy the called subroutine into `contains` section
2. Remove module dependencies
3. Pass any module variables as explicit arguments

### Unexecuted Kernels (Count=0)

For 181 kernels not executed in `test_real`:
- May need different simulation configuration to trigger
- Or mark as "data unavailable" and skip for now
- Can create synthetic test data as placeholder

---

## Build System Integration

### Compiling Dump Module

Add to main Makefile:
```makefile
DUMP_MOD = Kernel_benchmark/dump_kernel_data.o

$(DUMP_MOD): Kernel_benchmark/dump_kernel_data.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Add $(DUMP_MOD) to link step
solver.exe: $(OBJECTS) $(DUMP_MOD)
	$(FC) $(FFLAGS) -o $@ $^
```

### Conditional Compilation (Optional)

Use preprocessor to enable/disable dumping:
```fortran
#ifdef KERNEL_DUMP
  if (dump_call_count == DUMP_TARGET_CALL) then
    call dump_init(...)
    ...
  end if
#endif
```

Compile with: `$(FC) -DKERNEL_DUMP ...`

---

## Quality Checklist

For each completed benchmark:

- [ ] Benchmark compiles without errors
- [ ] Benchmark runs without crashes
- [ ] Validation passes (error count = 0)
- [ ] Timing is reasonable (compared to profiling data)
- [ ] All input arrays are read from files
- [ ] All output arrays are validated
- [ ] README.md updated with actual results
- [ ] Makefile works correctly

---

## Reference Files

| File | Description |
|------|-------------|
| `test_real/omp_profile.txt` | Profiling data (call counts, timing) |
| `meta_info_summary.md` | GPU difficulty analysis |
| `Kernel_benchmark/dump_kernel_data.f90` | Dump module |
| `Kernel_benchmark/102_diver3d_s_diver3d/` | Complete example |

---

## Expected Output

Successful benchmark run:
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

*Document created: 2026-01-02*
*Related: meta_info_instruction.md, meta_info_instruction2.md*
