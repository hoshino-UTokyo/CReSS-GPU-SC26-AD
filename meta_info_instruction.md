# GPU Porting Plan for MPI + OpenMP Weather Simulation (Claude Code Instructions)

## Ultimate Goal
The final goal is to port an existing **MPI + OpenMP** weather simulation codebase to run efficiently on **GPUs**.

## Project-Wide Phases (High-Level Roadmap)
1. **meta_info phase**
   - For every OpenMP `parallel` region (especially `parallel for` loops), perform **GPU-oriented loop analysis**.
   - Record findings directly in the source files as structured comments (see “Annotation Format” below).

2. **benchmark phase**
   - Extract every OpenMP `parallel` region into a **standalone benchmark** that can run independently.
   - Ensure each benchmark is minimal, reproducible, and representative of the original kernel behavior.

3. **benchmark-kernel GPU porting phase**
   - Port each benchmark kernel to GPU using **OpenACC** (initially).
   - Validate correctness and measure performance.

4. **simulation integration phase**
   - Integrate GPU-ported kernels back into the full simulation.
   - Enable end-to-end execution on GPU(s), including MPI decomposition + GPU kernels.

---

## What to Do in This Session (meta_info Phase)
You are currently working on the **meta_info phase**.

### Objective
For **every OpenMP parallel loop / region** in the codebase, determine the **GPU-porting difficulty** and identify blockers.

### For each OpenMP `parallel` region, check at minimum:
- **Thread-ID / OpenMP runtime dependencies**
  - e.g., usage of `omp_get_thread_num()` / `omp_get_thread_id()` or any logic depending on thread IDs.
- **Function calls inside the parallel region**
  - Identify whether calls are trivially inlinable / pure, or involve side effects / I/O.
- **Writes to global / shared state**
  - Writing to global variables, module-level variables, static state, or other shared memory.
- **Synchronization constructs**
  - Presence of `atomic`, `critical`, `ordered`, locks, reductions, barriers, etc.
- **Any other GPU blockers or hazards**
  - Indirect addressing, pointer aliasing, non-contiguous memory access, complex control flow, data dependencies, race risks,
    temporaries requiring privatization, implicit shared variables, etc.

### Output required
For each OpenMP `parallel` region, write an annotation **immediately above** that region in the source file, using exactly the format below.

---

## Annotation Format (MUST FOLLOW)
Insert this block right above the OpenMP parallel region:

!@llm start meta_info ----------------------------------------------------
! <your analysis here>
!@llm end meta_info ------------------------------------------------------


### Guidance for the analysis text
Make it concise but actionable. Include:
- **Location**: file name + function/subroutine name (if not obvious from context)
- **Loop summary**: what the loop does (1–2 lines)
- **GPU difficulty**: e.g., Easy / Medium / Hard
- **Reasons / blockers**: bullet list of specific issues found (thread-id use, atomics, global writes, function calls, etc.)
- **Suggested refactor** (if applicable): e.g., “convert global write to reduction”, “remove critical by privatization”, “inline function X”, etc.

Example:

!@llm start meta_info ----------------------------------------------------
! Location: foo.f90 :: subroutine update_state
! Summary : Updates prognostic variables over (i,j,k) domain.
! GPU diff: Medium
! Findings:
! - No omp_get_thread_* usage.
! - Calls calc_flux(...): appears pure, candidate for inlining / device routine.
! - Writes to module global diag_sum (race risk) -> convert to reduction or per-block buffer.
! - Uses atomic update on tend(i,j) -> consider privatization + reduction.
! Next:
! - Refactor diag_sum into local reduction.
! - Evaluate calc_flux for OpenACC routine eligibility.
!@llm end meta_info ------------------------------------------------------



---

## Workflow Constraints
- **Do not change program behavior** in this phase.
- Only add `meta_info` annotations and, if absolutely needed, tiny non-functional clarifications (e.g., comments).
- Make sure you cover **all** OpenMP parallel regions in the repository (search for `#pragma omp parallel`, `!$omp parallel`, `!$omp parallel do`, etc.).

---

## Git Requirements
You are on the **`meta_info` branch**.
After completing annotations:
1. Review changes (`git status`, `git diff`)
2. Commit with a clear message, e.g.:
   - `meta_info: annotate OpenMP regions for GPU porting`
