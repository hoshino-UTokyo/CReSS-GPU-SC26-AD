# CReSS GPU Porting Instructions

---

## IMPORTANT: Before Starting Any Work

**You MUST read these before starting work:**

1. **Latest progress report** - Understand current state and next steps
2. **Scripts reference** - Avoid recreating existing scripts

```bash
# Find and read the latest progress report
ls -t Claude/instructions/progress/progress_*.md | head -1

# Check existing scripts before writing new ones
cat Claude/scripts/README.md
```

**Checklist:**
1. Read the latest `progress/progress_YYYY-MM-DD.md` file
2. Review "Next Steps" section to understand what to work on
3. Check "In Progress" and "Blockers" sections for context
4. Check `Claude/scripts/README.md` for existing scripts

**At the end of each session:**

1. Create/update a progress report:
   ```bash
   cp Claude/instructions/progress/TEMPLATE.md Claude/instructions/progress/progress_YYYY-MM-DD.md
   ```

2. If you created new scripts, update `Claude/scripts/README.md`

See [Progress Report Guidelines](#progress-report-guidelines) below.

---

## Overview

This document provides instructions for porting the CReSS weather simulation codebase from **MPI + OpenMP** to **MPI + OpenACC (GPU)**.

## Target Architecture

- **Programming Model**: OpenACC for GPU offloading
- **Data Management**: NVIDIA Unified Memory (no explicit data directives)
- **Compiler**: nvfortran with `-acc -gpu=managed`

---

## Project Phases

The GPU porting process consists of 5 phases. Each phase has detailed instructions in the `tasks/` directory.

| Phase | Task | Description | Instruction |
|-------|------|-------------|-------------|
| 1 | Meta-Info Analysis | Analyze all OpenMP regions for GPU compatibility | [01_meta_info_analysis.md](tasks/01_meta_info_analysis.md) |
| 2 | Runtime Profiling | Add profiling to measure execution time and call counts | [02_runtime_profiling.md](tasks/02_runtime_profiling.md) |
| 3 | CPU Benchmark | Extract standalone benchmarks for each kernel | [03_cpu_benchmark.md](tasks/03_cpu_benchmark.md) |
| 4 | GPU Benchmark | Convert CPU benchmarks to GPU using OpenACC | [04_gpu_benchmark.md](tasks/04_gpu_benchmark.md) |
| 5 | Simulation Integration | Integrate GPU kernels back into main simulation | [05_simulation_integration.md](tasks/05_simulation_integration.md) |
| 6 | Kernel Optimization | nsysルーフライン解析に基づくカーネル最適化 | [06_kernel_optimization.md](tasks/06_kernel_optimization.md) |

---

## Quick Start

### Which phase are you working on?

1. **Analyzing OpenMP code for GPU compatibility?**
   - Read [01_meta_info_analysis.md](tasks/01_meta_info_analysis.md)
   - Output: Annotations added to source files

2. **Adding profiling to measure performance?**
   - Read [02_runtime_profiling.md](tasks/02_runtime_profiling.md)
   - Output: `omp_profile.txt` with timing data

3. **Creating standalone CPU benchmarks?**
   - Read [03_cpu_benchmark.md](tasks/03_cpu_benchmark.md)
   - Output: `Kernel_benchmark/<id>_<name>/`

4. **Converting to GPU with OpenACC?**
   - Read [04_gpu_benchmark.md](tasks/04_gpu_benchmark.md)
   - Output: `Kernel_benchmark_gpu/<id>_<name>/`

5. **Integrating GPU kernels into main simulation?**
   - Read [05_simulation_integration.md](tasks/05_simulation_integration.md)
   - Output: GPU-enabled `Src/*.f90` files

6. **Optimizing GPU kernel performance?**
   - Read [06_kernel_optimization.md](tasks/06_kernel_optimization.md)
   - Output: `Kernel_benchmark_gpu_opt/<id>_<name>/`

---

## Key Statistics

| Metric | Value |
|--------|-------|
| Total OpenMP parallel sections | 387 |
| Source files with OpenMP | 338 |
| GPU Difficulty: Easy | 236 (61%) |
| GPU Difficulty: Medium | 121 (31%) |
| GPU Difficulty: Hard | 30 (8%) |

For detailed analysis results, see [references/meta_info_summary.md](references/meta_info_summary.md).

---

## Directory Structure

```
CReSS/
├── Src/                          # Original source code (with meta_info annotations)
├── Kernel_benchmark/             # CPU benchmark programs
│   ├── dump_kernel_data.f90      # Shared dump module
│   ├── *.sh, *.py                # Benchmark scripts
│   └── <id>_<kernel_name>/       # Individual kernel benchmarks
├── Kernel_benchmark_gpu/         # GPU benchmark programs (OpenACC) — Phase 4 ベースライン
│   ├── Makefile.common           # Shared compiler settings
│   ├── build_all.sh, run_all.sh  # Build/run scripts
│   └── <id>_<kernel_name>/       # Individual GPU benchmarks
├── Kernel_benchmark_gpu_opt/     # 最適化版GPU benchmark — Phase 6
│   ├── roofline_analysis.sh      # nsys/ncuルーフライン解析
│   ├── optimization_summary.csv  # 全カーネル最適化結果
│   └── <id>_<kernel_name>/       # 最適化版 + OPTIMIZATION.md
├── test_real/                    # Real case simulation
│   ├── *.sh                      # Simulation run scripts
│   └── omp_profile.txt           # Profiling results
└── Claude/
    ├── instructions/             # Task documentation
    │   ├── main.md               # This file (READ FIRST)
    │   ├── tasks/                # Detailed task instructions
    │   ├── references/           # Analysis summaries
    │   └── progress/             # Session progress reports
    └── scripts/                  # Scripts reference
        └── README.md             # Catalog of all scripts
```

---

## GPU Porting Principles

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

### 2. Rely on Unified Memory

- No `!$acc data` directives needed
- Compile with `-gpu=managed`
- CUDA runtime handles data movement automatically

### 3. Explicit Loop Independence

- Always add `!$acc loop independent` to parallelizable loops
- Use `!$acc loop reduction(+:var)` for reduction loops

---

## Reference Documents

| Document | Purpose |
|----------|---------|
| [Meta Info Summary](references/meta_info_summary.md) | GPU porting difficulty analysis results |
| [Scripts Reference](../scripts/README.md) | Catalog of existing scripts |
| [omp_profile.txt](../../test_real/omp_profile.txt) | Runtime profiling data |

---

## Progress Report Guidelines

### Location

```
Claude/instructions/progress/
├── TEMPLATE.md              # Template for new reports
├── progress_2026-01-29.md   # Example: report for Jan 29
├── progress_2026-01-30.md   # Example: report for Jan 30
└── ...
```

### When to Create/Update

- **Start of session**: Read the latest progress report
- **End of session**: Create or update progress report for today

### Required Sections

| Section | Purpose |
|---------|---------|
| **Current Phase** | Which of the 4 phases are you in? |
| **Summary** | Brief description of session accomplishments |
| **Completed Tasks** | What was finished |
| **In Progress** | What is partially done |
| **Blockers / Issues** | Problems encountered |
| **Next Steps** | Priority list for next session (MOST IMPORTANT) |
| **Statistics** | Quantitative progress (kernels done, etc.) |

### Naming Convention

```
progress_YYYY-MM-DD.md
```

Example: `progress_2026-01-29.md`

### Tips

- Be specific in "Next Steps" - include file names, kernel IDs, etc.
- If multiple sessions in one day, append time: `progress_2026-01-29_14-30.md`
- Reference specific kernel directories: `Kernel_benchmark_gpu/041_bruntv_s_bruntv/`
- Include error messages in "Blockers" for debugging

---

*Last updated: 2026-03-30*
