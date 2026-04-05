# Progress Report: 2026-02-27 Session 2

## Phase 5: Simulation Integration - Verification & Additional Bug Fix

### Summary

After fixing phvuvw.f90, phvbcuvw.f90, lbcw.f90 in Session 1, proper comparison with CPU reference logs revealed that:
1. The **OLD GPU code (before fixes) already matched CPU** within ~0.0015% relative error at step 102
2. Our **fixes INTRODUCED ~0.125% divergence** at step 102
3. An additional staggered-grid array index bug was found and fixed, but had no effect

### Comparison Methodology

Three logs compared:
- **CPU reference**: `log.solver.20260226223516.txt` — compiled with `compile.conf_cpu_backup2` (no `-DUSE_GPU`, pure CPU), 361 steps
- **GPU_bug (bsearch_F1)**: `log.solver.bsearch_F1.txt` — old GPU code, 102 steps
- **GPU_fix1**: `log.solver.txt` (overwritten by fix2) — Session 1 fixes applied
- **GPU_fix2**: `log.solver.gpu_fix2.txt` — Session 2 additional index fix applied

Key configs:
- `compile.conf` = `compile.conf_bsearch_F1` (identical, verified via diff)
- `compile.conf_cpu_backup2`: FFLAGS = `-fast -mp -Mpreprocess -Mbyteswapio -mcmodel=medium` (no `-acc`, no `-DUSE_GPU`)

### TKE Comparison Results

| Step | CPU | GPU_bug | GPU_fix1 | GPU_fix2 | Bug vs CPU | Fix vs CPU |
|------|-----|---------|----------|----------|------------|------------|
| 4 | 0.56034401E-02 | 0.56034448E-02 | 0.56034448E-02 | 0.56034448E-02 | 0.0008% | 0.0008% |
| 40 | 0.74205095E+00 | 0.74204385E+00 | 0.74204379E+00 | 0.74204379E+00 | 0.001% | 0.001% |
| 83 | 0.39792508E-01 | 0.39793108E-01 | 0.39842319E-01 | 0.39842319E-01 | 0.0015% | **0.125%** |
| 102 | 0.32914385E-01 | 0.32914888E-01 | 0.32955598E-01 | 0.32955598E-01 | 0.0015% | **0.125%** |

- Steps 1-36: All nearly identical
- Steps 37-82: GPU_fix starts diverging slightly
- Steps 83-102: GPU_fix diverges ~0.125% from CPU, while GPU_bug stays ~0.0015%
- GPU_fix2 = GPU_fix1 (index fix had NO effect)

### Additional Bug Found and Fixed (No Effect)

Staggered-grid array index bugs in phvuvw.f90 GPU block mfcopt clamping:

| Location | GPU (wrong) | CPU (correct) | Fixed? | Active? |
|----------|-------------|---------------|--------|---------|
| v east (line 619) | `rmf8v(nim1,j,2)` | `rmf8v(nim2,j,2)` | Yes | No (exbvar(2:2)='-', not 'x') |
| w east (line 805) | `rmf(nim1,j,2)` | `rmf(nim2,j,2)` | Yes | Yes but no-op (phase speed < map factor) |
| v north (line 711) | `rmf8v(i,njm2,2)` | `rmf8v(i,njm1,2)` | Yes | No (exbvar(2:2)='-', not 'x') |

Root cause: v/w grids use different boundary indices than u-grid due to staggering.
These bugs are **real** (would matter for bc=4/5/6 or different exbvar configs) but have no effect on current config because:
- v not active in phvuvw (exbvar(2:2)='-')
- w east clamping is a no-op when constant phase speed < map factor

### Analysis: Why Fixes Introduced Divergence

The key mystery: OLD GPU code (no exbvar guards, no clamping, no dtsdb scaling) matched CPU, but NEW GPU code (with correct guards, clamping, scaling matching CPU logic) diverges.

**Phase speed data flow:**
```
phvuvw.f90 (called first) → writes phase speed arrays
phvbcuvw.f90 (called second) → overwrites phase speed arrays
```

**With exbvar='--x---xx':**
- phvuvw: only w active (exbvar(3:3)='x')
- phvbcuvw: u,v active (exbvar(1:1)='-', exbvar(2:2)='-'), w skipped

**OLD GPU code:**
1. phvuvw GPU: ALL components get constant (no guards) → u,v,w all set
2. phvbcuvw GPU: ALL components overwritten (no guards) → u,v,w all overwritten
3. Final: ALL phase speeds = phvbcuvw constant (unscaled by dtsdb)

**NEW GPU code:**
1. phvuvw GPU: only w gets constant + clamp + scale(×0.05)
2. phvbcuvw GPU: only u,v get constant + clamp + scale(×0.05)
3. Final: w = phvuvw scaled, u/v = phvbcuvw scaled

**CPU code:**
1. phvuvw CPU: only w gets constant + clamp + scale(×0.05)
2. phvbcuvw CPU: only u,v get constant + clamp + scale(×0.05)
3. Final: w = phvuvw scaled, u/v = phvbcuvw scaled

NEW GPU and CPU should produce **identical** final phase speeds. OLD GPU has **different** (unscaled) phase speeds.
Yet OLD GPU ≈ CPU and NEW GPU diverges 0.125%.

**Possible explanations:**
1. The divergence comes from a DIFFERENT GPU kernel (not 233/235/190), and OLD GPU's wrong phase speeds happened to cancel it (coincidence)
2. phvbcuvw (233) OLD GPU overwrites w too → the OLD code accidentally made w constant (matching some default), while NEW code correctly computes w from phvuvw → this correct computation has slightly different floating-point characteristics
3. The phvbcuvw GPU code (our new version) has a subtle bug not caught by the comparison — possibly in the complex phase speed calculation (bc=4/5/6 paths are not needed for bc=7, but maybe the exbvar guard logic is slightly different)

### Recommended Next Steps

1. **Revert phvuvw/phvbcuvw/lbcw to old GPU code** temporarily and verify that the 0.125% divergence disappears → this confirms which fix introduced the problem
2. **Isolate individual fixes**: Enable one fix at a time to identify which specific change causes the divergence
3. **Compare with all-GPU config**: The real test is whether the original 38% TKE divergence (with ALL kernels enabled) is fixed. The current bsearch_F1 config has 88 disabled kernels, and the observed 0.125% may be irrelevant to the main problem
4. **Try enabling ALL kernels** with current fixes to see if the 38% issue improves
5. **Fix remaining lower-priority bugs**: turbuvw.f90, strsten.f90, timeflt.f90

### Files Modified This Session

1. `Src/phvuvw.f90` — Fixed 3 staggered-grid array indices (nim1→nim2, njm2→njm1)

### Files Modified in Session 1 (same day)

1. `Src/phvuvw.f90` — Complete GPU block rewrite + loop bound fixes
2. `Src/phvbcuvw.f90` — Complete GPU block rewrite
3. `Src/lbcw.f90` — Loop bound fixes (8 occurrences)

### Key Log Files

| File | Description |
|------|-------------|
| `log.solver.20260226223516.txt` | CPU reference (all-CPU, 361 steps) |
| `log.solver.bsearch_F1.txt` | GPU_bug (old GPU code, 102 steps) |
| `log.solver.gpu_fix2.txt` | GPU after all fixes (this session, ~303+ steps) |
| `compile.conf` | Current build config (= bsearch_F1) |
| `compile.conf_cpu_backup2` | Pure CPU build config |
| `compile.conf_bsearch_F1` | Backup of bsearch_F1 config |
