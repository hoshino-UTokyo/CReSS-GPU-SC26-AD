# Phase 4: GPU Benchmark (OpenACC)

[Back to Main](../main.md)

## Objective

Convert CPU benchmarks from `Kernel_benchmark/` to GPU benchmarks using OpenACC in `Kernel_benchmark_gpu/`.

---

## Key Design Principles

### 1. Use OpenACC `kernels` Directive

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

**Why `kernels` instead of `parallel`?**
- Compiler analyzes and optimizes loops automatically
- Safer approach for complex loop nests
- Easier to port: just wrap existing loops

### 2. Unified Memory (No Data Directives)

**CRITICAL: Do NOT use `!$acc data` directives.**

Compile with `-gpu=managed` to enable NVIDIA Unified Memory.

### 3. Explicit Loop Independence

**CRITICAL: Always add `!$acc loop independent` to parallelizable loops.**

The compiler may be overly conservative. Explicit `independent` ensures parallelization.

---

## Conversion Process

### Step 0: Review CPU Benchmark README

**Read `Kernel_benchmark/<kernel_dir>/README.md` first!**

Check:
- **GPU Difficulty**: Prioritize Easy kernels
- **Findings**: Potential issues (reductions, subroutine calls)
- **Next**: Recommended porting approach

### Step 1: Copy CPU Benchmark

```bash
cp -r Kernel_benchmark/<kernel_dir> Kernel_benchmark_gpu/<kernel_dir>
```

### Step 2: Convert OpenMP to OpenACC

**Before (OpenMP):**
```fortran
!$omp parallel do private(i,j,k) schedule(runtime)
do k = 2, nk-2
  do j = 2, nj-2
    do i = 2, ni-2
      div(i,j,k) = (u(i+1,j,k) - u(i,j,k)) * dxiv
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
      div(i,j,k) = (u(i+1,j,k) - u(i,j,k)) * dxiv
    end do
  end do
end do
!$acc end kernels
```

### Step 3: Handle Special Cases

#### Multiple Loop Nests

Wrap each with separate `!$acc kernels`:

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

#### Reductions

Use `!$acc loop reduction(...)` (NOT `independent`):

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

#### Subroutine Calls

Add `!$acc routine seq` to called subroutines:

```fortran
!$acc kernels
!$acc loop independent
do k = 1, nk
  call compute_value(input(i,j,k), output(i,j,k))
end do
!$acc end kernels

contains

  !$acc routine seq
  subroutine compute_value(in_val, out_val)
    ...
  end subroutine
```

### Step 4: Update Makefile

Use the common Makefile include:

```makefile
include ../Makefile.common
```

`Kernel_benchmark_gpu/Makefile.common`:
```makefile
FC = nvfortran
FFLAGS = -O3 -acc -gpu=managed -mp -Mbyteswapio -Minfo=accel

BENCHMARK = kernel_benchmark

all: $(BENCHMARK)

$(BENCHMARK): kernel_benchmark.f90
	$(FC) $(FFLAGS) -o $@ $<

run: $(BENCHMARK)
	./$(BENCHMARK)

clean:
	rm -f $(BENCHMARK) *.o *.mod
```

### Step 5: Update Timing Code

Add `!$acc wait` for accurate timing:

```fortran
!$acc wait
t_start = omp_get_wtime()

call kernel_xxx(...)

!$acc wait
t_end = omp_get_wtime()
```

---

## Common Conversion Patterns

| Pattern | OpenMP | OpenACC |
|---------|--------|---------|
| Simple loop | `!$omp parallel do private(i,j,k)` | `!$acc kernels` + `!$acc loop independent` |
| Reduction | `!$omp parallel do reduction(+:sum)` | `!$acc loop reduction(+:sum)` |
| Multiple loops | `!$omp parallel` with multiple `!$omp do` | Separate `!$acc kernels` blocks |
| Subroutine call | (implicit) | Add `!$acc routine seq` |

---

## Directory Structure

```
Kernel_benchmark_gpu/
├── Makefile.common               # Shared compiler settings
├── build_all.sh                  # Build all benchmarks
├── run_all.sh                    # Run all benchmarks
│
├── 001_inidef_s_inidef/
│   ├── kernel_benchmark.f90      # GPU benchmark (OpenACC)
│   ├── Makefile                  # Includes ../Makefile.common
│   └── data -> ../../Kernel_benchmark/001_.../data  # Symlink
│
└── ... (kernel directories)
```

---

## Conversion Checklist

- [ ] Read `Kernel_benchmark/<kernel_dir>/README.md` first
- [ ] Copy directory to `Kernel_benchmark_gpu/`
- [ ] Create symlink for data directory
- [ ] Replace `!$omp parallel` with `!$acc kernels`
- [ ] Add `!$acc loop independent` to each parallelizable loop
- [ ] Remove OpenMP private/shared/schedule clauses
- [ ] Add `!$acc routine seq` to called subroutines
- [ ] Handle reductions with `!$acc loop reduction(...)`
- [ ] Create Makefile with `include ../Makefile.common`
- [ ] Add `!$acc wait` around timing calls
- [ ] Keep `use omp_lib` for timing
- [ ] Build and verify compilation
- [ ] Run and verify validation passes
- [ ] Check `-Minfo=accel` output
- [ ] Compare GPU time with CPU time (calculate speedup)

---

## Compiler Settings

### NVIDIA nvfortran (Recommended)

```bash
FC = nvfortran
FFLAGS = -O3 -acc -gpu=managed -mp -Minfo=accel

# Architecture-specific:
# -gpu=cc80    # A100
# -gpu=cc90    # H100
```

### Cray Compiler

```bash
FC = ftn
FFLAGS = -O3 -hacc -hmanaged -homp
```

---

## Debugging Tips

### Check Compiler Feedback

```bash
# Compile with -Minfo=accel shows what the compiler does:
nvfortran -Minfo=accel kernel_benchmark.f90
```

Output:
```
    123, Generating Tesla code
        125, !$acc loop gang, vector(128)
```

### Check for Data Movement

```bash
export NVCOMPILER_ACC_NOTIFY=1
./kernel_benchmark
```

### Common Errors

| Error | Solution |
|-------|----------|
| "Unsupported nested reduction" | Restructure or use atomics |
| "Accelerator region ignored" | Check loop bounds, function calls |
| "Call cannot be parallelized" | Add `!$acc routine seq` |

---

## Performance Notes

1. **Warmup is important** - JIT compilation on first execution
2. **Synchronization overhead** - `!$acc wait` has overhead
3. **Small kernels may not benefit** - Need enough parallelism
4. **Memory-bound** - Most weather kernels are memory-bound

---

## Expected Output

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

## Related Documents

- [Phase 3: CPU Benchmark](03_cpu_benchmark.md) - Source for conversion
- [Meta Info Summary](../references/meta_info_summary.md) - GPU difficulty info

---

*Original: gpu_kernel_instruction.md*
