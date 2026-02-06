# Runtime Profiling for OpenMP Parallel Sections (Claude Code Instructions)

## Purpose
Add runtime profiling to all OpenMP parallel sections to collect:
- **Execution count**: How many times each section is called
- **Loop length**: Total number of iterations per call
- **Execution time**: Wall-clock time spent in each section

This information helps identify hotspots and prioritize GPU porting efforts.

---

## Implementation Overview

### 1. Profiling Module (`Src/comprofile.f90`)

A Fortran module that provides:

```fortran
module m_comprofile
  ! Key functions:
  ! - profile_register(filename, subrname, description) -> section_id
  ! - profile_start(section_id)
  ! - profile_stop(section_id, loop_length)
  ! - profile_finalize()  ! Output results
end module
```

#### Module Features
- Supports up to 500 profiled sections
- Thread-safe using OpenMP atomics
- Uses `omp_get_wtime()` for high-resolution timing
- Outputs results to stdout and `omp_profile.txt`

---

### 2. Instrumentation Pattern

For each OpenMP parallel section, add the following:

#### Variable Declarations (in subroutine declaration section)
```fortran
! Profiling variables
integer, save :: prof_id1 = -1   ! Section ID (one per OMP section)
integer(8) :: loop_len           ! Loop length
```

#### Before `!$omp parallel`
```fortran
! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('filename.f90', 'subroutine_name', &
   & 'description')
end if
loop_len = int((k_end-k_start+1),8) * int((j_end-j_start+1),8) * int((i_end-i_start+1),8)
call profile_start(prof_id1)

!$omp parallel ...
```

#### After `!$omp end parallel`
```fortran
!$omp end parallel

call profile_stop(prof_id1, loop_len)
```

---

### 3. Module Reference

All files with OpenMP parallel sections must include:
```fortran
use m_comprofile
```

This should be added in the module reference section at the top of each file.

---

### 4. Program Finalization

In the main program (`solver.f90`), call `profile_finalize()` before program exit:

```fortran
! Output profiling results.
call profile_finalize()

! Finalize the MPI processes.
call endmpi(0)
```

---

## Output Format

The profiling output looks like:

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

========================================================================
```

---

## Build Instructions

### Compilation Order
`comprofile.f90` must be compiled before other source files that use it:

```bash
# Compile profiling module first
$FC $FFLAGS -c Src/comprofile.f90 -o Src/comprofile.o

# Then compile other sources
$FC $FFLAGS -c Src/other_file.f90 -o Src/other_file.o
```

### Required Flags
- OpenMP support: `-fopenmp` (gfortran) or `-qopenmp` (Intel)
- The module uses `omp_lib` for timing functions

### Linking
Include `comprofile.o` in the link step:
```bash
$FC $FFLAGS -o solver.exe Src/*.o -fopenmp
```

---

## Automation Scripts

Two Python scripts were created to automate instrumentation:

### `add_use_comprofile.py`
Adds `use m_comprofile` to all files with OpenMP parallel sections.

```bash
python3 add_use_comprofile.py
```

### `add_profiling_calls.py`
Adds profiling variable declarations and `profile_start`/`profile_stop` calls around each OpenMP parallel section.

```bash
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

### Loop Length Calculation
Loop length is calculated as the product of all loop bounds:
```fortran
loop_len = int((k_end-k_start+1),8) * int((j_end-j_start+1),8) * int((i_end-i_start+1),8)
```

For complex or conditional loops, a simplified estimate may be used.

### MPI Considerations
- Each MPI rank outputs its own profiling data
- Consider aggregating results across ranks for distributed runs
- The output file `omp_profile.txt` will be written by each rank (may overwrite)

### Disabling Profiling
To disable profiling without removing code:
1. Create a stub module that defines empty subroutines
2. Or use preprocessor directives to conditionally compile profiling code

---

## Future Enhancements

1. **MPI-aware output**: Aggregate results across MPI ranks
2. **Hierarchical profiling**: Track caller-callee relationships
3. **Memory usage**: Track memory allocation in each section
4. **GPU timing**: Extend to OpenACC kernel timing when ported
