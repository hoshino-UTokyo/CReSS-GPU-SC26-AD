# Meta Info Annotation Summary

## Overview

This document summarizes the GPU porting analysis annotations added to the CReSS codebase during the meta_info phase.

| Metric | Value |
|--------|-------|
| Annotated files | 338 |
| Total annotations | 387 |
| Lines added | 6,430 |
| Branch | `meta_info` |

## GPU Porting Difficulty Distribution

| Difficulty | Count | Percentage |
|------------|-------|------------|
| Easy | 236 | 61.0% |
| Medium | 121 | 31.3% |
| Hard | 30 | 7.8% |

### Difficulty Breakdown

```
Easy   [####################################          ] 236 (61.0%)
Medium [######################                        ] 121 (31.3%)
Hard   [#####                                         ]  30 ( 7.8%)
```

## Files with Multiple OpenMP Regions

| Annotations per file | File count |
|---------------------|------------|
| 1 | 305 |
| 2 | 26 |
| 3 | 4 |
| 4 | 1 |
| 6 | 1 |
| 8 | 1 |

## Hard Difficulty Files (30 files)

The "Hard" classification is based on algorithmic characteristics (vertical data dependencies,
sequential solvers, complex branching), **not** on OpenMP parallelism-level obstacles.
No structural porting blockers (atomic/critical, thread-ID dependencies, function calls
inside parallel regions, system function calls, etc.) were found in any file.

| File | Hard Factor (Algorithmic Characteristic) |
|------|------------------------------------------|
| `bulksfc.f90` | Complex physics with deep conditional branching |
| `coalbw.f90` | Bin microphysics coalescence |
| `defomten.f90` | Omega-theta computation |
| `depsitbw.f90` | Bin microphysics ice crystal growth |
| `exbcpt.f90` | External boundary conditions (multi-stage logic) |
| `exbcq.f90` | External boundary conditions (multi-stage logic) |
| `exbcss.f90` | External boundary conditions (multi-stage logic) |
| `exbcv.f90` | External boundary conditions (multi-stage logic) |
| `gaussel.f90` | Gaussian elimination (k-direction data dependency) |
| `get1d.f90` | 1D profile extraction (uses reduction) |
| `getkref.f90` | Reference level search (uses reduction) |
| `gseidel.f90` | Gauss-Seidel iteration (sequential dependency) |
| `inidisbw.f90` | Bin distribution initialization |
| `newblk_noevap.f90` | Bulk microphysics (complex conditional branching) |
| `phvbcs.f90` | Vertical physics boundary (4 directions x multi-stage) |
| `phvbcuvw.f90` | Vertical physics boundary (4 directions x multi-stage) |
| `phvs.f90` | Vertical physics scalar |
| `phvuvw.f90` | Vertical physics velocity |
| `rbcpt.f90` | Radiative boundary condition |
| `rbcq.f90` | Radiative boundary condition |
| `rbcqv.f90` | Radiative boundary condition |
| `rbcs0.f90` | Radiative boundary condition |
| `rbcs.f90` | Radiative boundary condition |
| `rbcw.f90` | Radiative boundary condition |
| `remapbw.f90` | Bin remapping |
| `setbin.f90` | Bin setup |
| `setblk.f90` | Block setup |
| `steps.f90` | Scalar time integration |
| `swadjst.f90` | SW adjustment |
| `vspdmp.f90` | Vertical sponge damping (uses reduction) |

## Files with OpenMP Reductions (25 files)

Files using `reduction()` clauses. These map directly to `!$acc loop reduction()`
in OpenACC and are **not** a fundamental porting obstacle.

- `adjstuv.f90`, `chkfile.f90`, `chkitr.f90`, `chkmoist.f90`, `chkmxn.f90`
- `cpondsfc.f90`, `fallblk.f90`, `fallbw.f90`, `fallqr.f90`, `get1d.f90`
- `getarea.f90`, `getkref.f90`, `getmxn.f90`, `hint2d.f90`, `hint3d.f90`
- `hintlnd.f90`, `newsindx.f90`, `outmxn.f90`, `paractl.f90`, `phycood.f90`
- `rdgrp.f90`, `set1d.f90`, `undefice.f90`, `undefsst.f90`, `vspdmp.f90`

## OpenACC Porting Obstacle Survey Results

A comprehensive survey of all OpenMP parallel regions found **no structural obstacles**
to OpenACC conversion:

| Check Item | Result |
|------------|--------|
| `omp_get_thread_num` / thread-ID dependency | **0 occurrences** — not used |
| `!$omp atomic` / `critical` / `ordered` | **0 occurrences** — not used |
| Function calls inside parallel regions (`call`) | **0 occurrences** — all calls are outside parallel regions |
| System/runtime functions (`malloc`, `free`, `system`, `getenv`, etc.) | **0 occurrences** — not used |
| Writes to global/module variables | **None** — module access is read-only |
| `reduction()` clauses | **25 files** — directly supported by `!$acc loop reduction()` |

In summary, all OpenMP parallel regions in CReSS are straightforward loop-level parallelism
with no structural obstacles for OpenMP-to-OpenACC conversion.
The "Hard" difficulty classification is based solely on algorithmic characteristics
(vertical dependencies, sequential solvers, complex branching).

## Common Patterns Found

### Easy (236 files)
- Pure stencil operations
- No function calls inside parallel regions
- No synchronization primitives
- Embarrassingly parallel grid loops
- Only intrinsic functions (max, min, sqrt, exp, log)

### Medium (121 files)
- Conditional branching inside loops
- Read-only access to module constants (`comphy`, `commath`, etc.)
- Multiple sequential parallel regions
- Helper function calls outside parallel regions

### Hard (30 files)
- k-direction data dependencies (vertical solvers: gaussel, gseidel)
- Complex conditional branching and multi-stage logic (boundary conditions, microphysics)
- Reduction clauses (directly supported by OpenACC)
- Complex control flow

## Annotation Format

Each annotation follows this structure:

```fortran
!@llm start meta_info ----------------------------------------------------
! Location: [filename] :: [subroutine name]
! Summary : [1-2 line description of what the region does]
! GPU diff: [Easy|Medium|Hard]
! Findings:
!   - [bullet points about thread-ID, function calls, synchronization]
! Next:
!   - [recommended GPU porting approach]
!@llm end meta_info ------------------------------------------------------
```

## Recommended Porting Priority

### Phase 1: Easy files (236 files)
- Direct conversion with `!$acc kernels` + `collapse`
- Minimal code changes required

### Phase 2: Medium files (121 files)
- May require kernel splitting for conditional branches
- Data region management considerations

### Phase 3: Hard files (30 files)
- Primarily algorithmic adaptation (no structural blockers)
- gaussel/gseidel: k-direction sequential dependency — parallelize over i,j with `!$acc kernels` (k remains sequential)
- Boundary conditions / microphysics: complex branching — handled by compiler via `!$acc kernels`
- Reductions: directly supported by `!$acc loop reduction()`

## Statistics by Category

| Category | Files | Description |
|----------|-------|-------------|
| Advection (`adv*`) | 6 | Advection operators |
| Adjustment (`adjst*`) | 10 | Field adjustment routines |
| Boundary (`bc*`, `rbc*`, `lbc*`, `exbc*`) | 25 | Boundary conditions |
| Turbulence (`turb*`) | 5 | Turbulence schemes |
| Surface (`sfc*`) | 8 | Surface physics |
| Microphysics | 40+ | Cloud/precipitation physics |
| Dynamics | 30+ | Core dynamics |
| I/O and Setup | 20+ | Initialization and output |

---

*Generated: 2025-12-28, Updated: 2026-04-04*
*Branch: meta_info*
*Commit: 4ddf2dc*
