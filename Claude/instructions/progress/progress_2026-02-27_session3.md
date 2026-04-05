# Progress Report: 2026-02-27 Session 3

## Phase 5: Simulation Integration - Compiler Flag Effect Isolation

### Summary

Extended the GPU fix2 log analysis beyond 102 steps, revealing catastrophic divergence at step 122+. Through systematic kernel categorization, safe-kernel testing, and all-disabled baseline testing, discovered that **all observed divergence comes from the `-acc -gpu=managed` compiler flags**, not from GPU kernels themselves. Established a proper GPU-flagged baseline for future kernel comparison.

### Extended GPU fix2 Analysis (303 steps)

Three divergence phases identified when comparing gpu_fix2 (66 GPU kernels) vs CPU reference:

| Phase | Steps | TKE Relative Error | Pattern |
|-------|-------|--------------------|---------|
| Phase 1 | 1-55 | ±0.004% | Small oscillating |
| Phase 2 | 56-121 | +0.1252% | Stable positive bias |
| Phase 3 | 122+ | -47.16% | Catastrophic drop |

Step 122 = 2 steps after dump at step 120 (dmpitv=600s). Kernel 327 (timeflt, Asselin filter) is DISABLED in bsearch_F1.

### Kernel Safety Categorization

Analyzed all 66 enabled kernels in bsearch_F1 config and categorized:

| Category | Count | Criteria |
|----------|-------|----------|
| SAFE | 28 | Simple data-parallel, no branching, straightforward index mapping |
| MEDIUM | 25 | Moderate complexity, conditional branches, multiple arrays |
| RISKY | 13 | Known bugs, very complex, or large kernels |

**28 SAFE kernels**: 181, 183, 190, 213, 214, 231, 237, 247, 264, 265, 266, 278, 285, 288, 293, 294, 302, 303, 305, 315, 321, 332, 333, 353, 358, 360, 361, 362

### Three-Way Comparison Test

Built and ran three configurations, all compared against CPU reference (`log.solver.20260226223516.txt`):

| Config | GPU Kernels | Compiler Flags | TKE at step 83 | Catastrophe |
|--------|-------------|----------------|-----------------|-------------|
| CPU reference | 0 | `-fast -mp` (no `-acc`) | baseline | none |
| safe_only | 28 | `-fast -mp -acc -gpu=managed` | +0.125% | -28% at step 153 |
| all-disabled | 0 | `-fast -mp -acc -gpu=managed` | +0.125% | -28% at step 153 |
| gpu_fix2 (bsearch_F1) | 66 | `-fast -mp -acc -gpu=managed` | +0.125% | -47% at step 122 |

### Critical Finding

**safe_only and all-disabled produce IDENTICAL results** (byte-for-byte same TKE values at every step).

This proves:
1. The 28 SAFE GPU kernels introduce **zero additional error** beyond compiler flag effects
2. **ALL divergence** (±0.004%, +0.1252%, -28%) comes from the `-acc -gpu=managed` compiler flags alone
3. These flags change floating-point behavior for ALL code, including CPU-executed paths (optimization differences, FP contraction, etc.)

### OpenACC Directive Leak Check

Verified that no `!$acc` directives exist outside `#ifdef USE_GPU` guards:
- Preprocessed all source files with `cpp -E` using all DISABLE flags
- Found **0 surviving `!$acc` directives** in preprocessed output
- All divergence is confirmed to be from compiler flag effects, not stray OpenACC directives

### Remaining Kernel Issue

The gpu_fix2 run (66 kernels) shows -47% at step 122, while all-disabled/safe_only shows -28% at step 153. The earlier and more severe catastrophe means **some of the 38 additional kernels (MEDIUM/RISKY) DO cause extra divergence** beyond compiler flag effects. These need investigation via binary search.

### Files Created This Session

| File | Description |
|------|-------------|
| `compile.conf_safe_only` | Config with 28 SAFE GPU kernels enabled (128 DISABLE flags) |
| `compile.conf_allcpu_accflag` | Config with ALL 156 kernels disabled but `-acc -gpu=managed` flags |
| `Src/solver_gpu_all_disabled.exe` | Binary compiled with all-disabled config |
| `test_real/log.solver.gpu_all_disabled.txt` | Reference log from all-disabled run (303 steps) |
| `compare_tke.sh`, `compare_tke_safe.sh`, `compare_tke_full.sh`, `compare_tke_detail.sh` | TKE comparison scripts |

### Key Reference Logs

| File | Config | Steps | Purpose |
|------|--------|-------|---------|
| `log.solver.20260226223516.txt` | Pure CPU (no `-acc`) | 361 | Absolute CPU reference |
| `log.solver.gpu_all_disabled.txt` | All-disabled + `-acc -gpu=managed` | 303 | GPU-flagged baseline (compiler effect only) |
| `log.solver.gpu_fix2.txt` | bsearch_F1 (66 kernels) | 303 | GPU with fixes applied |

### Recommended Next Steps

1. **Investigate compiler flag divergence**: Try `-acc` without `-gpu=managed`, or use `-Kieee` to enforce strict FP
2. **Use all-disabled log as proper baseline**: Future kernel comparisons should use `log.solver.gpu_all_disabled.txt` instead of pure CPU reference, to isolate kernel effects from compiler flag effects
3. **Binary search among 38 MEDIUM/RISKY kernels**: Identify which kernels cause the earlier catastrophe (step 122 vs step 153)
4. **Consider `-Mnofma` or `-Kieee`**: These flags may reduce FP divergence from compiler optimization differences
