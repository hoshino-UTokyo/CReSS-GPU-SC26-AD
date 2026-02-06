# Phase 5: Simulation Integration

[Back to Main](../main.md)

## Objective

Integrate GPU-ported kernels from `Kernel_benchmark_gpu/` back into the main CReSS simulation codebase (`Src/`). Enable end-to-end GPU execution with MPI decomposition.

---

## Prerequisites

Before starting integration:

- [ ] All target GPU benchmarks pass validation (error count = 0)
- [ ] GPU benchmarks show acceptable performance
- [ ] Integration target list is prepared

---

## Integration Strategy

### Batch Integration with Binary Search Debugging

**Rationale**: Full simulation takes ~2 hours per run. Testing each kernel individually is not practical. Since benchmarks already validate correctness, apply all kernels at once and use binary search if issues arise.

```
1. Apply ALL GPU kernels at once
2. Run full simulation test
3. If PASS → Done
4. If FAIL → Binary search to identify problematic kernel(s)
```

### Binary Search Process (When Issues Occur)

```
All kernels (N) fail
    ↓
Enable first half (N/2) only → Test
    ↓
If PASS: Problem in second half
If FAIL: Problem in first half
    ↓
Repeat until single kernel identified
```

**Example with 100 kernels:**
```
Run 1: All 100 kernels      → FAIL
Run 2: Kernels 1-50         → PASS  (problem in 51-100)
Run 3: Kernels 1-50 + 51-75 → FAIL  (problem in 51-75)
Run 4: Kernels 1-50 + 51-62 → PASS  (problem in 63-75)
Run 5: Kernels 1-50 + 51-62 + 63-69 → FAIL (problem in 63-69)
Run 6: Kernels 1-50 + 51-62 + 63-66 → PASS (problem in 67-69)
Run 7: Kernels 1-50 + 51-62 + 63-66 + 67-68 → PASS (problem is kernel 69)
```

With 100 kernels, only ~7 test runs needed (log2(100) ≈ 7) instead of 100.

### Preprocessor-Based Kernel Control

Each kernel is controlled by preprocessor macros:

```fortran
#if defined(USE_GPU) && !defined(DISABLE_GPU_091)
  ! OpenACC version
#else
  ! OpenMP version (original code preserved)
#endif
```

**Macro definitions:**

| Macro | Effect |
|-------|--------|
| `USE_GPU` | Enable GPU code paths (master switch) |
| `DISABLE_GPU_091` | Force kernel 091 to use OpenMP |
| `DISABLE_GPU_102` | Force kernel 102 to use OpenMP |
| ... | (one per kernel) |

**Compile-time control:**

```bash
# All GPU
-DUSE_GPU

# All CPU (omit USE_GPU)
(no flags)

# GPU except kernels 091, 102
-DUSE_GPU -DDISABLE_GPU_091 -DDISABLE_GPU_102
```

This design keeps original OpenMP code intact while allowing flexible GPU/CPU switching per kernel.

---

## Step-by-Step Integration Process

### Step 1: Prepare Integration List

Create a list of all kernels to integrate:

```bash
# List all validated GPU benchmarks
ls Kernel_benchmark_gpu/*/kernel_benchmark.f90 | wc -l

# Create integration list
cat > integration_list.txt << 'EOF'
# Kernel ID | Source File | Status
091_diver3d | Src/diver3d.f90 | pending
102_pgrad | Src/pgrad.f90 | pending
...
EOF
```

### Step 2: Add OpenACC with Preprocessor Guards

**IMPORTANT: Keep existing OpenMP code intact. Add OpenACC as an alternative controlled by preprocessor.**

For each kernel in the list, wrap with preprocessor and add OpenACC version:

**Before (Original OpenMP - keep this):**
```fortran
!$omp parallel do private(i,j,k) schedule(runtime)
do k = 2, nk-2
  do j = 2, nj-2
    do i = 2, ni-2
      div(i,j,k) = (u(i+1,j,k) - u(i,j,k)) * dxiv &
                 + (v(i,j+1,k) - v(i,j,k)) * dyiv
    end do
  end do
end do
!$omp end parallel do
```

**After (OpenMP preserved + OpenACC added with preprocessor):**
```fortran
#if defined(USE_GPU) && !defined(DISABLE_GPU_091)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
!$acc kernels
!$acc loop independent
do k = 2, nk-2
  !$acc loop independent
  do j = 2, nj-2
    !$acc loop independent
    do i = 2, ni-2
      div(i,j,k) = (u(i+1,j,k) - u(i,j,k)) * dxiv &
                 + (v(i,j+1,k) - v(i,j,k)) * dyiv
    end do
  end do
end do
!$acc end kernels
#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel do private(i,j,k) schedule(runtime)
do k = 2, nk-2
  do j = 2, nj-2
    do i = 2, ni-2
      div(i,j,k) = (u(i+1,j,k) - u(i,j,k)) * dxiv &
                 + (v(i,j+1,k) - v(i,j,k)) * dyiv
    end do
  end do
end do
!$omp end parallel do
#endif
```

**Design benefits:**
- Original OpenMP code remains intact (safe rollback)
- GPU/CPU switching at compile time
- Individual kernel disable for binary search debugging
- Single source file maintains both versions

### Step 3: Handle Subroutine Calls

If kernels call subroutines, add `!$acc routine` directive:

```fortran
!$acc routine seq
subroutine helper_function(...)
  ...
end subroutine
```

**Location options:**
1. Add directive in the original subroutine file
2. Create a separate interface block if modifying original is not possible

### Step 4: Update Build System

#### Makefile Changes

Add OpenACC flags and preprocessor support:

```makefile
# Compiler
FC = nvfortran

# Base flags
FFLAGS = -O3 -mp -Mbyteswapio

# GPU flags (OpenACC + Unified Memory)
FFLAGS_GPU = -acc -gpu=managed -Minfo=accel

# Preprocessor flags (passed via command line or here)
# USE_GPU: Enable GPU code paths
# DISABLE_GPU_XXX: Disable specific kernel XXX for binary search
CPPFLAGS = -DUSE_GPU

# Combined flags
FFLAGS_ALL = $(FFLAGS) $(FFLAGS_GPU) $(CPPFLAGS) $(EXTRA_FLAGS)

# Build rule
%.o: %.f90
    $(FC) $(FFLAGS_ALL) -c $< -o $@
```

#### Compile Options

```bash
# All GPU enabled (normal operation)
make CPPFLAGS="-DUSE_GPU"

# All CPU / OpenMP only (baseline, no GPU)
make CPPFLAGS=""

# Disable specific kernels for binary search
make CPPFLAGS="-DUSE_GPU" EXTRA_FLAGS="-DDISABLE_GPU_091 -DDISABLE_GPU_102"

# Disable range of kernels (e.g., 051-100)
DISABLE_FLAGS=$(for i in $(seq 51 100); do echo -n "-DDISABLE_GPU_$(printf '%03d' $i) "; done)
make CPPFLAGS="-DUSE_GPU" EXTRA_FLAGS="$DISABLE_FLAGS"
```

### Step 5: Build and Verify Compilation

```bash
make clean
make 2>&1 | tee build.log

# Check for OpenACC messages
grep -i "Generating\|acc\|kernel" build.log
```

Expected output:
```
diver3d.f90:
    123, Generating Tesla code
        125, !$acc loop gang, vector(128)
```

### Step 6: Run Full Simulation Test

```bash
cd test_real
./run_test.sh   # Takes ~2 hours

# Compare output with reference
diff output.dat reference/output.dat
```

### Step 7: If Test Fails - Binary Search

If simulation produces wrong results:

```bash
# Create binary search script
cat > binary_search_kernels.sh << 'EOF'
#!/bin/bash
# Disable kernels in second half (51-100)
DISABLE_FLAGS=""
for i in $(seq 51 100); do
    DISABLE_FLAGS="$DISABLE_FLAGS -DDISABLE_GPU_$(printf '%03d' $i)"
done
make clean
make EXTRA_FLAGS="$DISABLE_FLAGS"
./run_test.sh
EOF
```

Track binary search progress:
```markdown
| Run | Enabled Kernels | Result | Conclusion |
|-----|-----------------|--------|------------|
| 1   | All (1-100)     | FAIL   | Problem exists |
| 2   | 1-50            | PASS   | Problem in 51-100 |
| 3   | 1-75            | FAIL   | Problem in 51-75 |
| ...
```

### Step 8: Fix and Re-test

Once problematic kernel identified:

1. Review GPU benchmark for that kernel
2. Check for differences between benchmark and integrated code
3. Fix the issue
4. Re-run binary search to find any additional problems

### Step 9: Performance Verification (After All Tests Pass)

```bash
# Run CPU version (baseline)
export ACC_DEVICE_TYPE=host
time ./solver.exe

# Run GPU version
export ACC_DEVICE_TYPE=nvidia
time ./solver.exe

# Calculate speedup
```

---

## Data Management Considerations

### Unified Memory Approach (Recommended)

With `-gpu=managed`, no explicit data directives needed:

```fortran
! Data automatically migrates between CPU and GPU
!$acc kernels
do k = 1, nk
  ...
end do
!$acc end kernels
! Data available on CPU after kernel
```

### Explicit Data Regions (Advanced)

For performance optimization, consider explicit data regions:

```fortran
!$acc data copyin(input1, input2) copyout(output)

!$acc kernels
do k = 1, nk
  ...
end do
!$acc end kernels

!$acc kernels
do k = 1, nk
  ...
end do
!$acc end kernels

!$acc end data
```

**When to use explicit data:**
- Multiple kernels operating on same arrays
- Large arrays that should stay on GPU
- MPI communication boundaries (need explicit sync)

---

## MPI + GPU Considerations

### GPU-Aware MPI

If using GPU-aware MPI (e.g., CUDA-aware OpenMPI):

```fortran
! Data can be sent directly from GPU memory
!$acc host_data use_device(sendbuf, recvbuf)
call MPI_Sendrecv(sendbuf, ..., recvbuf, ...)
!$acc end host_data
```

### Standard MPI (Non GPU-Aware)

Data must be on host for MPI:

```fortran
! Ensure data is on host before MPI
!$acc update host(sendbuf)
call MPI_Sendrecv(sendbuf, ..., recvbuf, ...)
!$acc update device(recvbuf)
```

### Halo Exchange Pattern

```fortran
! 1. Compute interior (GPU)
!$acc kernels
do k = 2, nk-1
  do j = 2, nj-1
    do i = 2, ni-1
      ...
    end do
  end do
end do
!$acc end kernels

! 2. Update host for MPI
!$acc update host(array)

! 3. MPI halo exchange (CPU)
call exchange_halo(array, ...)

! 4. Update device with new halo
!$acc update device(array)

! 5. Continue computation (GPU)
!$acc kernels
...
!$acc end kernels
```

---

## Debugging Integration Issues

### Common Problems

| Problem | Cause | Solution |
|---------|-------|----------|
| Wrong results | Race condition | Check `independent` clause is correct |
| Slow performance | Excessive data transfer | Add explicit data regions |
| MPI hangs | GPU/CPU data mismatch | Add `!$acc update` before MPI |
| Compilation error | Missing routine directive | Add `!$acc routine seq` |

### Debug Environment Variables

```bash
# Show data transfer
export NVCOMPILER_ACC_NOTIFY=1

# Show kernel launches
export NVCOMPILER_ACC_NOTIFY=2

# Detailed debug info
export NVCOMPILER_ACC_NOTIFY=3

# Synchronous execution (easier debugging)
export NVCOMPILER_ACC_SYNCHRONOUS=1
```

### Validation Script

```bash
#!/bin/bash
# validate_integration.sh

# Run CPU version
export ACC_DEVICE_TYPE=host
./solver.exe
mv output.dat output_cpu.dat

# Run GPU version
export ACC_DEVICE_TYPE=nvidia
./solver.exe
mv output.dat output_gpu.dat

# Compare
python3 compare_outputs.py output_cpu.dat output_gpu.dat
```

---

## Integration Checklist

### Pre-Integration (Per Kernel)

- [ ] GPU benchmark validation passed (error count = 0)
- [ ] Benchmark README reviewed for any special notes

### Batch Integration

- [ ] All target kernels have passing GPU benchmarks
- [ ] Integration list created (`integration_list.txt`)
- [ ] All OpenMP directives replaced with OpenACC
- [ ] All `!$acc routine seq` added to called subroutines
- [ ] Preprocessor guards added for binary search capability
- [ ] Makefile updated with `-acc -gpu=managed`

### Validation

- [ ] Compilation succeeds with accelerator info messages
- [ ] Full simulation test completed (~2 hours)
- [ ] Output matches reference data
- [ ] If failed: binary search completed, problems identified and fixed

### Final

- [ ] All kernels integrated and validated
- [ ] Performance measured (GPU vs CPU speedup)
- [ ] Progress report updated with final status

---

## Tracking Integration Status

Create/update `integration_status.md` in progress folder:

```markdown
# Integration Status

## Summary
- Total kernels to integrate: X
- Integrated: Y
- Validation: PASS/FAIL
- Binary search runs: Z

## Kernel Status

| Kernel ID | Name | Benchmark OK | Integrated | Issue Found | Notes |
|-----------|------|--------------|------------|-------------|-------|
| 091 | diver3d | Yes | Yes | No | |
| 102 | pgrad | Yes | Yes | No | |
| 115 | advuvw | Yes | Yes | Yes | Fixed race condition |
| 120 | turbs | Yes | Yes | No | |

## Binary Search Log (if applicable)

| Run | Date | Kernels Enabled | Result | Next Action |
|-----|------|-----------------|--------|-------------|
| 1 | 01-29 | All (1-100) | FAIL | Test 1-50 |
| 2 | 01-29 | 1-50 | PASS | Test 51-75 |
| ...
```

---

## Rollback Procedure

### Quick Rollback (Preprocessor)

```bash
# Compile without GPU (all kernels use OpenMP)
make clean
make USE_GPU=0

# Or disable specific kernels
make EXTRA_FLAGS="-DDISABLE_GPU_091 -DDISABLE_GPU_102"
```

### Full Rollback (Git)

```bash
# Revert all source changes
git checkout HEAD -- Src/

# Or revert specific file
git checkout HEAD -- Src/diver3d.f90
```

---

## Performance Optimization (After Integration)

Once basic integration works:

1. **Profile GPU execution**
   ```bash
   nsys profile ./solver.exe
   ncu --target-processes all ./solver.exe
   ```

2. **Identify bottlenecks**
   - Data transfer overhead
   - Kernel launch overhead
   - Memory access patterns

3. **Optimize data movement**
   - Add explicit data regions
   - Overlap computation and transfer
   - Use async clauses

4. **Optimize kernels**
   - Adjust gang/worker/vector
   - Consider loop tiling
   - Optimize memory access

---

## Related Documents

- [Phase 4: GPU Benchmark](04_gpu_benchmark.md) - Source of GPU code
- [Meta Info Summary](../references/meta_info_summary.md) - Kernel dependencies
- [Progress Reports](../progress/) - Track integration progress

---

*Created: 2026-01-29*
