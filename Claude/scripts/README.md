# Scripts Reference

This document catalogs existing scripts in the project for easy reference across sessions.

---

## Profiling Scripts (Phase 2)

| Script | Location | Purpose |
|--------|----------|---------|
| `add_use_comprofile.py` | Root | Add `use m_comprofile` to source files |
| `add_profiling_calls.py` | Root | Add profiling instrumentation to OpenMP sections |
| `add_profiling.py` | Root | Profiling helper |
| `add_runtime_info.py` | Root | Add runtime info to source |

---

## Data Dump Scripts (Phase 3)

| Script | Location | Purpose |
|--------|----------|---------|
| `add_dump_code.py` | Root | Add dump instrumentation to source files |
| `extract_all_variables.py` | Root | Extract variable lists from OpenMP sections |
| `check_dump_coverage.py` | Root | Check which kernels have dump data |
| `fix_dump_fp_values.py` | Root | Fix floating point values in dump |

---

## CPU Benchmark Scripts (Phase 3)

Location: `Kernel_benchmark/`

| Script | Purpose |
|--------|---------|
| `build_all_benchmarks.sh` | Build all CPU benchmarks |
| `run_all_benchmarks.sh` | Run all CPU benchmarks |
| `check_data.sh` | Check data files exist |
| `check_dump_coverage.py` | Verify dump coverage |
| `copy_dump_data.sh` | Copy dump data to benchmark dirs |
| `ln_dump_data.sh` | Create symlinks for dump data |
| `setup_data_symlinks.sh` | Setup data directory symlinks |
| `find_incomplete.sh` | Find benchmarks missing data |
| `extract_dump_info.sh` | Extract dump information |
| `update_readme_with_dump_info.sh` | Update README with dump info |
| `update_readme_dump.py` | Update README (Python version) |

---

## GPU Benchmark Scripts (Phase 4)

Location: `Kernel_benchmark_gpu/`

| Script | Purpose |
|--------|---------|
| `build_all.sh` | Build all GPU benchmarks |
| `run_all.sh` | Run all GPU benchmarks |

---

## Simulation Scripts

Location: `test_real/`

| Script | Purpose |
|--------|---------|
| `solver.sh` | Run solver |
| `MIYABI-solver.sh` | Run solver on MIYABI |
| `MIYABI-solver-dump.sh` | Run solver with dump enabled |
| `MIYABI-pripro.sh` | Run preprocessor on MIYABI |
| `MIYABI-unite.sh` | Run unite on MIYABI |
| `FLOW-*.sh` | Various FLOW scripts |

Location: `test_ideal/`

| Script | Purpose |
|--------|---------|
| `FLOW-unite.sh` | Unite for ideal test |
| `FLOW-small.sh` | Small ideal test |

---

## Utility Scripts

| Script | Location | Purpose |
|--------|----------|---------|
| `create_kernel_dirs.py` | Root | Create kernel benchmark directories |
| `fix_missing_use.py` | Src/ | Fix missing use statements |
| `ispack_download.sh` | Root | Download ISPACK library |

---

## Usage Examples

### Build and run all CPU benchmarks

```bash
cd Kernel_benchmark
./build_all_benchmarks.sh
./run_all_benchmarks.sh
```

### Build and run all GPU benchmarks

```bash
cd Kernel_benchmark_gpu
./build_all.sh
./run_all.sh
```

### Check dump data coverage

```bash
python3 check_dump_coverage.py
```

### Run simulation with dump

```bash
cd test_real
./MIYABI-solver-dump.sh
```

---

## Adding New Scripts

When creating new scripts:

1. Place in appropriate location (or `Claude/scripts/` if general-purpose)
2. Update this README with description
3. Add usage example if non-trivial

---

*Last updated: 2026-01-29*
