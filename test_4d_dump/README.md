# 4D Array Dump Segfault Analysis

## Problem Summary

Solver crashes with segmentation fault during 4D array dump in `steps.f90:591`.

**Error messages:**
```
Caught signal 11 (Segmentation fault: invalid permissions for mapped object)
Cgroup memsw limit exceeded
```

## Root Cause

**Cumulative I/O buffer memory exhaustion.**

Analysis of `test_real/log.solver.txt` shows:
- **753 3D arrays** were dumped before the crash
- Each 3D array: ~104M elements (~396 MB)
- 4D array dump (`qwtrf_in.bin`) crashes at ~648th dump cycle

The NVIDIA Fortran compiler's I/O library accumulates internal buffers during stream writes. Even with explicit `flush()` and `close()`, memory is not fully released. After ~650 large array writes, the cumulative memory usage exceeds the cgroup memory+swap limit.

## Reproduction

```bash
cd test_4d_dump
make
./test_cumulative_dump 700   # Will fail around dump 648
./test_cumulative_dump 100   # Will succeed
./test_4d_dump               # Single dump test - always succeeds
```

## Solution Options

### Option 1: Batch Dumping (Recommended)

Dump kernels in batches across multiple solver runs instead of all at once.

**Implementation:**
1. Add environment variable or config to select which kernels to dump:
```fortran
! At module level
integer, save :: dump_batch = 0  ! 0=no dump, 1=batch1, 2=batch2, ...

! At init
character(len=10) :: env_val
call get_environment_variable('DUMP_BATCH', env_val)
read(env_val, *, iostat=ios) dump_batch
```

2. Run solver multiple times with different DUMP_BATCH values:
```bash
DUMP_BATCH=1 ./solver.exe   # Dumps kernels 1-50
DUMP_BATCH=2 ./solver.exe   # Dumps kernels 51-100
...
```

### Option 2: Selective Dumping

Only dump the specific kernels you need for benchmarking.

**Implementation:**
Add a configuration file listing kernels to dump:
```
# dump_targets.conf
steps
diverpih
advp
```

### Option 3: Memory-Optimized Dump Module

Modify dump functions to minimize memory usage:

```fortran
subroutine dump_array_4d_safe(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
  ! Write element-by-element to avoid temp array creation
  integer :: i, j, k, l

  open(...)
  do l = l1, l2
    do k = k1, k2
      do j = j1, j2
        do i = i1, i2
          write(dump_unit) arr(i, j, k, l)
        end do
      end do
    end do
  end do
  close(...)
end subroutine
```

This is slower but avoids the ~3MB temp array created by `arr(i1:i2, j1:j2, k, l)`.

### Option 4: Periodic Sync

Force OS buffer flush periodically:

```fortran
! After every N dumps
if (mod(dump_count, 50) == 0) then
  call execute_command_line('sync', wait=.true.)
end if
```

## Recommended Action

**Use Option 1 (Batch Dumping)** as the primary solution:

1. Modify `dump_kernel_data.f90` to support batch selection
2. Create a script that runs solver multiple times with different batches
3. Each run dumps ~50-100 kernels to stay well under the memory limit

This requires minimal code changes and avoids the memory issue entirely.

## Files

- `test_4d_dump.f90` - Single 4D dump test (demonstrates basic functionality works)
- `test_cumulative_dump.f90` - Cumulative dump test (reproduces the crash)
- `Makefile` - Build configuration
