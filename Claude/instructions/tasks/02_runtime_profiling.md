# Phase 2: Runtime Profiling

[Back to Main](../main.md)

## Objective

Add runtime profiling to all OpenMP parallel sections to collect:
- **Execution count**: How many times each section is called
- **Loop length**: Total number of iterations per call
- **Execution time**: Wall-clock time spent in each section

This data helps identify hotspots and prioritize GPU porting efforts.

---

## Implementation Overview

### Profiling Module

Location: `Src/comprofile.f90`

```fortran
module m_comprofile
  ! Key functions:
  ! - profile_register(filename, subrname, description) -> section_id
  ! - profile_start(section_id)
  ! - profile_stop(section_id, loop_length)
  ! - profile_finalize()  ! Output results
end module
```

---

## Instrumentation Pattern

### Step 1: Add Module Reference

```fortran
use m_comprofile
```

### Step 2: Add Variables (in declaration section)

```fortran
! Profiling variables
integer, save :: prof_id1 = -1   ! Section ID (one per OMP section)
integer(8) :: loop_len           ! Loop length
```

### Step 3: Add Profiling Calls

```fortran
! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('filename.f90', 'subroutine_name', &
   & 'description')
end if
loop_len = int((k_end-k_start+1),8) * int((j_end-j_start+1),8) * int((i_end-i_start+1),8)
call profile_start(prof_id1)

!$omp parallel ...
...
!$omp end parallel

call profile_stop(prof_id1, loop_len)
```

### Step 4: Add Finalization (in main program)

```fortran
! In solver.f90, before program exit
call profile_finalize()
```

---

## Output Format

Results are written to `omp_profile.txt`:

```
========================================================================
OpenMP Parallel Section Profiling Results
========================================================================

ID    File                 Subroutine           Count        Avg Loops   Total Time(s)   Avg Time(ms)
------------------------------------------------------------------------------------------------------------------------
    1 gsmoow.f90           s_gsmoow             1000         1280000.0        1.234567        1.234567
    2 gsmoow.f90           s_gsmoow             1000         1000000.0        0.987654        0.987654
    3 steppe.f90           s_steppe             5000          640000.0        2.345678        0.469136
...
```

---

## Build Instructions

### Compilation Order

`comprofile.f90` must be compiled first:

```bash
# Compile profiling module first
$FC $FFLAGS -c Src/comprofile.f90 -o Src/comprofile.o

# Then compile other sources
$FC $FFLAGS -c Src/other_file.f90 -o Src/other_file.o
```

### Linking

```bash
$FC $FFLAGS -o solver.exe Src/*.o -fopenmp
```

---

## Automation Scripts

| Script | Purpose |
|--------|---------|
| `add_use_comprofile.py` | Adds `use m_comprofile` to all files |
| `add_profiling_calls.py` | Adds profiling calls around OpenMP sections |

```bash
python3 add_use_comprofile.py
python3 add_profiling_calls.py
```

---

## Current Status

| Item | Count |
|------|-------|
| Files with `use m_comprofile` | 298 |
| Files with profiling calls | 339 |
| Total OpenMP sections instrumented | ~387 |

---

## Notes

### MPI Considerations

- Each MPI rank outputs its own profiling data
- Output file may be overwritten by different ranks
- Consider aggregating results across ranks for distributed runs

### Disabling Profiling

1. Create a stub module with empty subroutines, or
2. Use preprocessor directives for conditional compilation

---

## Related Documents

- [Phase 1: Meta-Info Analysis](01_meta_info_analysis.md) - Run before this phase
- [Phase 3: CPU Benchmark](03_cpu_benchmark.md) - Uses profiling data

---

*Original: meta_info_instruction2.md*
