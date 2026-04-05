# Phase 1: Meta-Info Analysis

[Back to Main](../main.md)

## Objective

Analyze every OpenMP parallel region in the codebase and determine GPU porting difficulty. Add structured annotations directly to source files.

---

## What to Check for Each OpenMP Region

| Check Item | Description |
|------------|-------------|
| Thread-ID dependencies | Usage of `omp_get_thread_num()`, `omp_get_thread_id()` |
| Function calls | Calls inside parallel region (pure? side effects? I/O?) |
| Global/shared writes | Writes to module variables, global state |
| Synchronization | `atomic`, `critical`, `ordered`, locks, barriers |
| Reductions | `reduction(+:var)` clauses |
| System/runtime calls | `malloc`, `free`, `system`, `getenv`, `exit`, `abort` etc. inside parallel region |
| Other hazards | Indirect addressing, non-contiguous access, control flow |

---

## Annotation Format

Insert this block **immediately above** each OpenMP parallel region:

```fortran
!@llm start meta_info ----------------------------------------------------
! Location: [filename] :: [subroutine name]
! Summary : [1-2 line description]
! GPU diff: [Easy|Medium|Hard]
! Findings:
!   - [bullet points about issues found]
! Next:
!   - [recommended GPU porting approach]
!@llm end meta_info ------------------------------------------------------
```

### Example

```fortran
!@llm start meta_info ----------------------------------------------------
! Location: foo.f90 :: subroutine update_state
! Summary : Updates prognostic variables over (i,j,k) domain.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - Calls calc_flux(...): appears pure, candidate for inlining / device routine.
!   - Writes to module global diag_sum (race risk) -> convert to reduction or per-block buffer.
!   - Uses atomic update on tend(i,j) -> consider privatization + reduction.
! Next:
!   - Refactor diag_sum into local reduction.
!   - Evaluate calc_flux for OpenACC routine eligibility.
!@llm end meta_info ------------------------------------------------------

!$omp parallel do private(i,j,k) ...
```

---

## GPU Difficulty Classification

### Easy

- Pure stencil operations
- No function calls inside parallel regions
- No synchronization (atomic, critical, barrier)
- Embarrassingly parallel grid loops
- Only intrinsic functions (max, min, sqrt, exp, log)

**Porting approach**: Direct `!$acc kernels` with `loop independent`

### Medium

- Conditional branching inside loops
- Module constant access
- Multiple sequential parallel regions
- Calls to simple helper functions

**Porting approach**: May need function inlining or `!$acc routine seq`

### Hard

- Thread-ID dependent operations
- Reduction operations
- k-dependent data dependencies (vertical solver)
- Gauss elimination/Gauss-Seidel solvers
- Complex control flow with early exits

**Porting approach**: Algorithm redesign may be required

---

## Workflow

1. Search for all `!$omp parallel` regions in `Src/`
2. For each region, analyze using the checklist above
3. Add annotation block immediately above the region
4. **Do not change program behavior** - only add annotations
5. Generate a summary report (see below)

---

## Report Output

After annotating all regions, generate a summary report at `Claude/instructions/references/meta_info_summary.md` containing:

1. **Statistics**: Total number of OpenMP regions, breakdown by difficulty (Easy/Medium/Hard)
2. **Per-file table**: File name, subroutine, difficulty, key findings (1-line each)
3. **Hazard inventory**: List of all regions with system/runtime calls, reductions, synchronization, thread-ID dependencies
4. **Recommended porting order**: Suggested order based on difficulty and dependency relationships

This report serves as the planning input for subsequent phases.

---

## Git Workflow

```bash
# Ensure you're on the meta_info branch
git checkout meta_info

# After completing annotations
git add Src/*.f90
git commit -m "meta_info: annotate OpenMP regions for GPU porting"
```

---

## Related Documents

- [Phase 2: Runtime Profiling](02_runtime_profiling.md) - Add after this phase
- [Meta Info Summary](../references/meta_info_summary.md) - Analysis results

---

*Original: meta_info_instruction.md*
