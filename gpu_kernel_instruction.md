# GPU Kernel Benchmark Guide (OpenACC Version)

## Purpose

Create GPU-accelerated versions of each kernel benchmark using OpenACC. This guide describes how to convert CPU benchmarks in `Kernel_benchmark/` to GPU benchmarks in `Kernel_benchmark_gpu/`.

---

## Key Design Principles

### 1. OpenACC as the Programming Model

Use OpenACC directives for GPU acceleration instead of OpenMP.

**Rationale:**
- OpenACC provides portable GPU acceleration for Fortran
- Directive-based approach allows incremental porting
- Compiler handles GPU memory management and kernel launches

### 2. Unified Memory (No Explicit Data Directives)

**CRITICAL: Do NOT use `!$acc data` directives.**

Instead, rely on NVIDIA Unified Memory (managed memory) to handle CPU-GPU data transfers automatically.

**Rationale:**
- Simplifies code by eliminating explicit data movement
- Reduces porting effort and potential bugs
- Enables gradual migration without restructuring data flow
- Performance is sufficient for benchmarking purposes

**Compiler Flags Required:**
```bash
-gpu=managed          # Enable Unified Memory
-Minfo=accel          # Show accelerator information
```

### 3. Use `kernels` Directive Instead of `parallel`

**CRITICAL: Use `!$acc kernels` instead of `!$acc parallel`.**

```fortran
! CORRECT: Use kernels directive
!$acc kernels
do k = 1, nk
  do j = 1, nj
    do i = 1, ni
      output(i,j,k) = input(i,j,k) * factor
    end do
  end do
end do
!$acc end kernels

! AVOID: parallel directive with explicit loop
!$acc parallel loop collapse(3)
do k = 1, nk
  ...
```

**Rationale:**
- `kernels` directive lets the compiler analyze and optimize loops automatically
- Compiler determines parallelism, loop scheduling, and gang/worker/vector mapping
- More conservative but safer approach for complex loop nests
- Easier to port: just wrap existing loop with `!$acc kernels ... !$acc end kernels`

### 4. Explicit Loop Parallelization with `!$acc loop independent`

**CRITICAL: Do NOT leave everything to the compiler. Use `!$acc loop independent` to explicitly indicate parallelizable loops.**

```fortran
!$acc kernels
!$acc loop independent
do k = 1, nk
  !$acc loop independent
  do j = 1, nj
    !$acc loop independent
    do i = 1, ni
      output(i,j,k) = input(i,j,k) * factor
    end do
  end do
end do
!$acc end kernels
```

**Rationale:**
- Compiler may be overly conservative and fail to parallelize safe loops
- `independent` clause explicitly tells the compiler "this loop has no loop-carried dependencies"
- Ensures intended parallelization without relying on compiler analysis
- Makes parallelization intent clear in the code

**When to use `independent`:**
- Loops where each iteration is completely independent (no cross-iteration dependencies)
- Typical stencil operations, element-wise operations, array assignments

**When NOT to use `independent`:**
- Loops with reduction operations (use `!$acc loop reduction(...)` instead)
- Loops with actual dependencies between iterations
- Sequential loops that must not be parallelized

---

## Conversion Process: CPU to GPU Benchmark

### Step 0: Review CPU Benchmark README

**CRITICAL: Before starting conversion, read `Kernel_benchmark/<kernel_dir>/README.md`.**

The README contains valuable information for GPU porting:

```
Kernel_benchmark/<kernel_dir>/README.md
├── Source Location      - Original file and line number
├── GPU Difficulty       - Easy/Medium/Hard rating
├── Summary              - What the kernel does
├── Runtime Profile      - CPU execution time (baseline for comparison)
│   ├── Calls            - Number of calls during simulation
│   ├── Average Loop Length
│   ├── Total Time
│   └── Average Time per Call
└── Original Annotation  - Detailed analysis
    ├── Findings         - Key observations (dependencies, special cases)
    └── Next             - Recommended GPU porting approach
```

**Pay special attention to:**
- **GPU Difficulty**: Prioritize "Easy" kernels first
- **Findings**: Lists potential issues (e.g., `omp_get_thread_*` usage, subroutine calls, reductions)
- **Next**: Specific recommendations for OpenACC porting (e.g., "collapse(3)", "separate kernels for branches")
- **Average Time per Call**: Use as baseline to evaluate GPU speedup

**Example from README:**
```
## Analysis
- **GPU Difficulty**: Easy
- **Summary**: Forces hydrometeor mixing ratios to be non-negative

## Original Annotation
! Findings:
!   - No omp_get_thread_* usage.
!   - Pure max() operations, all GPU compatible.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i).
```

### Step 1: Copy CPU Benchmark Directory

```bash
cp -r Kernel_benchmark/<kernel_dir> Kernel_benchmark_gpu/<kernel_dir>
```

### Step 2: Convert OpenMP to OpenACC

Replace OpenMP directives with OpenACC `kernels` directive:

**Before (OpenMP):**
```fortran
!$omp parallel do private(i,j,k) schedule(runtime)
do k = 2, nk-2
  do j = 2, nj-2
    do i = 2, ni-2
      div(i,j,k) = (u(i+1,j,k) - u(i,j,k)) * dxiv &
                 + (v(i,j+1,k) - v(i,j,k)) * dyiv &
                 + (w(i,j,k+1) - w(i,j,k)) * dziv
    end do
  end do
end do
!$omp end parallel do
```

**After (OpenACC):**
```fortran
!$acc kernels
!$acc loop independent
do k = 2, nk-2
  !$acc loop independent
  do j = 2, nj-2
    !$acc loop independent
    do i = 2, ni-2
      div(i,j,k) = (u(i+1,j,k) - u(i,j,k)) * dxiv &
                 + (v(i,j+1,k) - v(i,j,k)) * dyiv &
                 + (w(i,j,k+1) - w(i,j,k)) * dziv
    end do
  end do
end do
!$acc end kernels
```

### Step 3: Handle Special Cases

#### 3.1 Multiple Loop Nests

Wrap each independent loop nest with its own `!$acc kernels` region:

```fortran
!$acc kernels
!$acc loop independent
do k = 1, nk
  !$acc loop independent
  do j = 1, nj
    !$acc loop independent
    do i = 1, ni
      tmp1(i,j,k) = ...
    end do
  end do
end do
!$acc end kernels

!$acc kernels
!$acc loop independent
do k = 1, nk
  !$acc loop independent
  do j = 1, nj
    !$acc loop independent
    do i = 1, ni
      output(i,j,k) = tmp1(i,j,k) + ...
    end do
  end do
end do
!$acc end kernels
```

#### 3.2 Reductions

Use `!$acc loop reduction` for reduction loops (do NOT use `independent` on reduction loops):

```fortran
sum_val = 0.0
!$acc kernels
!$acc loop reduction(+:sum_val)
do k = 1, nk
  !$acc loop reduction(+:sum_val)
  do j = 1, nj
    !$acc loop reduction(+:sum_val)
    do i = 1, ni
      sum_val = sum_val + array(i,j,k)
    end do
  end do
end do
!$acc end kernels
```

**Note:** For reductions, use `reduction` clause instead of `independent`. The reduction clause implies parallelization with proper reduction handling.

#### 3.3 Conditional Execution Inside Loops

Conditionals inside loops work naturally with `kernels`. Use `independent` as usual:

```fortran
!$acc kernels
!$acc loop independent
do k = 1, nk
  !$acc loop independent
  do j = 1, nj
    !$acc loop independent
    do i = 1, ni
      if (mask(i,j,k) > 0) then
        output(i,j,k) = input(i,j,k) * factor
      else
        output(i,j,k) = 0.0
      end if
    end do
  end do
end do
!$acc end kernels
```

#### 3.4 Subroutine Calls Inside Kernels

If the kernel calls subroutines, add `!$acc routine seq` to the called subroutine:

```fortran
!$acc kernels
!$acc loop independent
do k = 1, nk
  !$acc loop independent
  do j = 1, nj
    !$acc loop independent
    do i = 1, ni
      call compute_value(input(i,j,k), output(i,j,k))
    end do
  end do
end do
!$acc end kernels

contains

  !$acc routine seq
  subroutine compute_value(in_val, out_val)
    real, intent(in) :: in_val
    real, intent(out) :: out_val
    out_val = in_val * 2.0
  end subroutine
```

### Step 4: Update Makefile

**Use the common Makefile include** - each benchmark only needs a simple Makefile:

```makefile
# GPU Kernel Benchmark Makefile
# Include common settings from parent directory
include ../Makefile.common
```

The common settings are defined in `Kernel_benchmark_gpu/Makefile.common`:

```makefile
# Common Makefile settings for GPU Kernel Benchmarks

# Compiler settings for NVIDIA GPU with OpenACC
FC = nvfortran
FFLAGS = -O3 -acc -gpu=managed -mp -Mbyteswapio -Minfo=accel

# Note on flags:
#   -O3           : Optimization level 3
#   -acc          : Enable OpenACC
#   -gpu=managed  : Use NVIDIA Unified Memory (managed memory)
#   -mp           : Enable OpenMP (for omp_get_wtime)
#   -Mbyteswapio  : Enable byte-swapping for binary I/O (required for big-endian data files)
#   -Minfo=accel  : Show accelerator optimization information

# Target executable name
BENCHMARK = kernel_benchmark

# Default target
all: $(BENCHMARK)

# Build rule
$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

# Run target
run: $(BENCHMARK)
	./$(BENCHMARK)

# Clean target
clean:
	rm -f $(BENCHMARK) *.o *.mod

.PHONY: all run clean
```

### Step 5: Update Timing Code

Keep using `omp_get_wtime()` for timing. Add `!$acc wait` to ensure GPU operations complete before timing.

```fortran
use omp_lib
...
! Ensure GPU operations complete before timing
!$acc wait
t_start = omp_get_wtime()

call kernel_xxx(...)

! Ensure GPU operations complete
!$acc wait
t_end = omp_get_wtime()
```

**Note:** OpenMP and OpenACC can coexist. Keep `-mp` (or `-fopenmp`) flag in Makefile alongside OpenACC flags.

---

## GPU Benchmark Program Template

```fortran
!***********************************************************************
! GPU Kernel Benchmark: <kernel_name>
!***********************************************************************
! Source: Src/<filename>.f90
! Subroutine: <subroutine_name>
! GPU Port: OpenACC with Unified Memory
!***********************************************************************
program kernel_benchmark_gpu_<name>
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

  ! 6. Warmup iterations (includes JIT compilation)
  do iter = 1, warmup_iterations
    call kernel_<name>(...)
    !$acc wait
  end do

  ! 7. Benchmark iterations (timed)
  t_total = 0.0d0
  do iter = 1, num_iterations
    !$acc wait
    t_start = omp_get_wtime()

    call kernel_<name>(...)

    !$acc wait
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
  ! The GPU kernel (converted from OpenMP to OpenACC)
  !-----------------------------------------------------------------
  subroutine kernel_<name>(...)
    ! === OpenACC kernels region ===
    !$acc kernels
    !$acc loop independent
    do k = ...
      !$acc loop independent
      do j = ...
        !$acc loop independent
        do i = ...
          ! Kernel computation
        end do
      end do
    end do
    !$acc end kernels
  end subroutine

  !-----------------------------------------------------------------
  ! I/O utilities (same as CPU version)
  !-----------------------------------------------------------------
  subroutine read_config(...)
  subroutine read_parameters(...)
  subroutine read_array_2d(...)
  subroutine read_array_3d(...)
  subroutine validate_output(...)
  subroutine print_results(...)

end program
```

---

## Directory Structure

```
Kernel_benchmark_gpu/
├── Makefile.common                # Common compiler settings (shared by all benchmarks)
├── build_all.sh                   # Script to build all benchmarks
├── run_all.sh                     # Script to run all benchmarks
├── README.md                      # Index of GPU kernels
│
├── 001_inidef_s_inidef/
│   ├── kernel_benchmark.f90       # GPU benchmark (OpenACC)
│   ├── Makefile                   # Includes ../Makefile.common
│   ├── benchmark.conf             # Runtime configuration
│   └── data -> ../../Kernel_benchmark/001_.../data  # Symlink to CPU data
│
├── 091_diver3d_s_diver3d/
│   ├── kernel_benchmark.f90
│   ├── Makefile
│   └── data -> ...
│
└── ... (kernel directories)
```

**Notes:**
- Use symlinks to share data files with CPU benchmarks to save disk space.
- Each benchmark's Makefile simply includes `../Makefile.common` to use shared compiler settings.

---

## Conversion Checklist

For each kernel being converted to GPU:

- [ ] **Read `Kernel_benchmark/<kernel_dir>/README.md` first** (check GPU Difficulty, Findings, Next)
- [ ] Copy directory from `Kernel_benchmark/` to `Kernel_benchmark_gpu/`
- [ ] Create symlink for data directory
- [ ] Replace `!$omp parallel` with `!$acc kernels`
- [ ] Add `!$acc loop independent` to each parallelizable loop
- [ ] Remove OpenMP private/shared/schedule clauses (not needed for kernels)
- [ ] Add `!$acc routine seq` to any called subroutines
- [ ] Handle reductions with `!$acc loop reduction(...)` (not `independent`)
- [ ] Create Makefile that includes `../Makefile.common`
- [ ] Add `!$acc wait` before and after timing with `omp_get_wtime()`
- [ ] Keep `use omp_lib` (OpenMP and OpenACC coexist)
- [ ] Build and verify compilation succeeds
- [ ] Run and verify validation passes
- [ ] Check compiler output for `-Minfo=accel` messages
- [ ] Compare GPU time with CPU time from README (calculate speedup)
- [ ] Update GPU benchmark README with GPU execution time and speedup

---

## Common Conversion Patterns

### Pattern 1: Simple Triple-Nested Loop

**OpenMP:**
```fortran
!$omp parallel do private(i,j,k)
do k = 1, nk
  do j = 1, nj
    do i = 1, ni
      a(i,j,k) = b(i,j,k) + c(i,j,k)
    end do
  end do
end do
!$omp end parallel do
```

**OpenACC:**
```fortran
!$acc kernels
!$acc loop independent
do k = 1, nk
  !$acc loop independent
  do j = 1, nj
    !$acc loop independent
    do i = 1, ni
      a(i,j,k) = b(i,j,k) + c(i,j,k)
    end do
  end do
end do
!$acc end kernels
```

### Pattern 2: Loop with Reduction

**OpenMP:**
```fortran
sum = 0.0
!$omp parallel do reduction(+:sum) private(i,j,k)
do k = 1, nk
  do j = 1, nj
    do i = 1, ni
      sum = sum + a(i,j,k)**2
    end do
  end do
end do
!$omp end parallel do
```

**OpenACC:**
```fortran
sum = 0.0
!$acc kernels
!$acc loop reduction(+:sum)
do k = 1, nk
  !$acc loop reduction(+:sum)
  do j = 1, nj
    !$acc loop reduction(+:sum)
    do i = 1, ni
      sum = sum + a(i,j,k)**2
    end do
  end do
end do
!$acc end kernels
```

### Pattern 3: Multiple Independent Loops

**OpenMP:**
```fortran
!$omp parallel
!$omp do private(i,j,k)
do k = 1, nk
  ...
end do
!$omp end do

!$omp do private(i,j,k)
do k = 1, nk
  ...
end do
!$omp end do
!$omp end parallel
```

**OpenACC:**
```fortran
!$acc kernels
!$acc loop independent
do k = 1, nk
  ...
end do
!$acc end kernels

!$acc kernels
!$acc loop independent
do k = 1, nk
  ...
end do
!$acc end kernels
```

### Pattern 4: Loop with Subroutine Call

**OpenMP:**
```fortran
!$omp parallel do private(i,j,k)
do k = 1, nk
  do j = 1, nj
    do i = 1, ni
      call helper(a(i,j,k), b(i,j,k))
    end do
  end do
end do
!$omp end parallel do

contains
  subroutine helper(x, y)
    ...
  end subroutine
```

**OpenACC:**
```fortran
!$acc kernels
!$acc loop independent
do k = 1, nk
  !$acc loop independent
  do j = 1, nj
    !$acc loop independent
    do i = 1, ni
      call helper(a(i,j,k), b(i,j,k))
    end do
  end do
end do
!$acc end kernels

contains
  !$acc routine seq
  subroutine helper(x, y)
    ...
  end subroutine
```

---

## Compiler-Specific Notes

### NVIDIA nvfortran (Recommended)

```bash
FC = nvfortran
FFLAGS = -O3 -acc -gpu=managed -mp -Minfo=accel

# -mp                # Enable OpenMP (for omp_get_wtime)
# -gpu=cc80          # Target specific GPU architecture (A100)
# -gpu=cc90          # Target specific GPU architecture (H100)
# -gpu=fastmath      # Enable fast math operations
# -Mcuda             # Enable CUDA interoperability
```

### Cray Compiler (ftn)

```bash
FC = ftn
FFLAGS = -O3 -hacc -hmanaged -homp

# -homp              # Enable OpenMP (for omp_get_wtime)
# -hnoacc            # Disable OpenACC (for comparison)
```

---

## Debugging Tips

### 1. Check Compiler Feedback

Always compile with `-Minfo=accel` to see what the compiler is doing:

```
kernel_benchmark.f90:
    123, Generating Tesla code
        125, !$acc loop gang, vector(128) ! blockIdx.x threadIdx.x
        127, !$acc loop gang ! blockIdx.y
```

### 2. Verify Kernel Execution

Add runtime checks:

```fortran
print *, 'Before kernel'
!$acc kernels
...
!$acc end kernels
!$acc wait
print *, 'After kernel'
```

### 3. Check for Data Movement

With Unified Memory, data movement happens automatically, but you can verify:

```bash
# Set environment variable to see data transfers
export NVCOMPILER_ACC_NOTIFY=1
./kernel_benchmark_gpu
```

### 4. Common Errors

**Error: "Unsupported nested reduction"**
- Solution: Restructure the reduction or use atomic operations

**Error: "Accelerator region ignored"**
- Check that loops have computable bounds
- Ensure no function calls that can't be inlined

**Error: "Call to ... cannot be parallelized"**
- Add `!$acc routine seq` to the called function

---

## Performance Considerations

### 1. Warmup is Important

GPU kernel compilation (JIT) happens on first execution. Always run warmup iterations before timing.

### 2. Synchronization Overhead

`!$acc wait` has overhead. For timing accuracy, it's necessary, but minimize in production code.

### 3. Small Kernels May Not Benefit

Very small loop nests may not have enough parallelism to benefit from GPU execution. Focus on large kernels first.

### 4. Memory-Bound vs Compute-Bound

Most weather simulation kernels are memory-bound. GPU performance depends heavily on memory bandwidth utilization.

---

## Expected Output

Successful GPU benchmark run:

```
==================================================
 GPU Kernel Benchmark: diver3d
==================================================
 Grid size: ni=   902, nj=   902, nk=   128
 GPU Device: 0 (NVIDIA A100-SXM4-40GB)
==================================================
 Average time:     2.341 ms
 Total time:      23.410 ms
 Min time:         2.298 ms
 Max time:         2.456 ms
--------------------------------------------------
 Max relative error:   1.2345E-07
 Tolerance:            1.0000E-05
 Error count:                   0
 Validation: PASSED
==================================================
```

---

## Reference

| Item | Description |
|------|-------------|
| CPU Benchmarks | `Kernel_benchmark/` |
| GPU Benchmarks | `Kernel_benchmark_gpu/` |
| CPU Instruction | `benchmark_instruction.md` |
| OpenACC Specification | [OpenACC 3.3](https://www.openacc.org/specification) |
| NVIDIA nvfortran | [HPC SDK Documentation](https://docs.nvidia.com/hpc-sdk/) |

---

*Document created: 2026-01-21*
*Related: benchmark_instruction.md*
