# CReSS GPU port — SC26 AD/AE Artifact $A_1$

This repository is **Artifact $A_1$** of the SC26 AD/AE submission for the
paper on AI-assisted GPU porting of the CReSS regional weather code.
It contains the GPU-ported CReSS source tree with the full Phase 1–6
human-authored specification (`Claude/instructions/`), per-kernel CPU
and GPU benchmark directories, and the TC2214 typhoon evaluation
scenario.

The companion artifact $A_2$ (Phase 3 controlled experiments) is in a
separate repository:
<https://github.com/hoshino-UTokyo/CReSS-GPU-SC26-AD-phase3>.

## Paper's main contributions

- **$C_1$** Reformulation of GPU porting of a large-scale legacy HPC
  application as an iterative human-in-the-loop workflow under
  CLI-based AI agent session constraints, rather than a one-shot code
  generation problem. The workflow is exercised end-to-end on the
  250,000+ line CReSS weather simulation code, with all 162 active
  OpenMP parallel regions ported to GPU through OpenACC.

- **$C_2$** Identification of the dominant instability factors in
  AI-assisted GPU porting workflows from the perspective of costly
  re-execution and lack of global cost awareness: (i) local task
  execution without cost awareness, (ii) incomplete variable
  specification leading to invalid dump acquisition, and (iii)
  session-induced state loss. Quantitatively demonstrated in the
  Phase 3 controlled experiments ($A_2$).

- **$C_3$** A set of specification components (workflow ordering,
  per-kernel `variable_list`, I/O configuration, and recovery rule)
  that externalize execution context and enable stable workflow
  convergence across session boundaries. The role of each component is
  isolated by the ablation study in $A_2$.

| Artifact | Contributions supported | Related paper elements |
|----------|------------------------|------------------------|
| $A_1$ (this repo) | $C_1$, $C_3$ | Sec. III, Sec. IV-D, Fig. 1 |
| $A_2$ (phase3 repo) | $C_2$, $C_3$ | Sec. IV-B (Tab. II), Sec. IV-C (Tab. III) |

Both repositories are released under a paper-verification-only license
(see `LICENSE`); a future, broadly-licensed open-source release of the
GPU port is planned separately. A persistent DOI for each repository
will be issued at the AE stage.

## Relation to contributions

$A_1$ is the GPU-ported CReSS application produced by the
human-in-the-loop workflow defined in Sec. III. It substantiates
$C_1$ by showing that the proposed workflow yields a complete GPU port
of a 250,000+ line production weather code, and $C_3$ by exposing the
specification components (`Claude/instructions/`, including
`tasks/01_meta_info_analysis.md` through `tasks/06_kernel_optimization.md`,
plus the per-kernel `variable_list.txt` and `benchmark.conf`) as a
concrete realization of externalized execution context.

## Expected results

Building $A_1$ on the target hardware reproduces the end-to-end GPU
port. Running the TC2214 simulation for 360 time steps on a single
Miyabi-G node yields the validation file
`test_real/result/TC2214.mon.check.txt`, whose maximum and minimum
pressure values $pp_{\max}$ and $pp_{\min}$ satisfy the correctness
criterion of Sec. III-A,
$|a_1| < 10^{-4}$ and $|a_2| < 10^{-4}$, with

$$a_1 = (pp_{\max} - 1.140663\times10^3) / 1.140663\times10^3$$
$$a_2 = (pp_{\min} + 3.927225\times10^3) / 3.927225\times10^3$$

The resulting GPU run achieves a median per-step execution time of
~1.88 s versus ~9.51 s on 72 Grace CPU OpenMP threads
(speedup ~5.1×, Sec. IV-D). Profiling all 162 active kernels with
`nsys` and the roofline analysis classifies them as memory-bound, with
the dominant kernels sustaining 35–60% of peak HBM bandwidth. The
top-10 kernels account for ~65% of total GPU time, reproducing Fig. 1.

## Reproduction time

- **Setup** ~60 min (load modules, build ISPACK, build CReSS, stage
  GPV input).
- **Execution** ~15 min for the 360-step GPU TC2214 run; ~60 min for
  the corresponding 72-thread CPU run used for the speedup baseline;
  ~30 min for full `nsys` profiling of all kernels.
- **Analysis** ~15 min (run validation script on
  `TC2214.mon.check.txt`, regenerate kernel ranking from `nsys` output).

## Hardware

A single GPU node of the Miyabi-G partition at JCAHPC, equipped with
one NVIDIA GH200 Grace–Hopper Superchip (Grace CPU, 72 cores @ 3.0 GHz,
120 GB system memory at ~512 GB/s; NVIDIA Hopper H100 GPU with 96 GB
HBM at ~4,022 GB/s, NVLink-C2C 450 GB/s cache-coherent interconnect).
Any functionally equivalent node with one GH200 is sufficient. The
node is used as an interactive job (queue `short-g`) with a 2 h
walltime limit; this constraint is essential because the session-reset
behavior of the workflow originates from this limit. For the speedup
baseline reported in Sec. IV-D the same node is used with the GPU
disabled (CPU-only build). About 200 GB of scratch storage is required
(89 GB dump data plus binaries and logs).

## Software

- Operating system: Red Hat Enterprise Linux 9.4 (`aarch64`, kernel
  5.14.0-427).
- NVIDIA HPC SDK 25.9 (`nvfortran`, NCCL, HPC-X)
  <https://developer.nvidia.com/hpc-sdk>.
- CUDA Toolkit bundled with HPC SDK 25.9.
- MPI: HPC-X (or Open MPI 4.1+).
- Profiling: `nsys` 2024.x, `ncu` (bundled with HPC SDK).
- ISPACK 0.95 (<https://www.gfd-dennou.org/arch/ispack/>); the
  bundled `ispack_download.sh` fetches the exact version required.
- Node.js LTS and `@anthropic-ai/claude-code` (Claude Code CLI).
  Version pinned to Opus 4.6 for the controlled experiments. Provided
  by the wrapper script `claude.sh`.
- Python 3.9+ for the helper scripts under `Claude/scripts/` and the
  project root (`add_dump_code.py`, `check_dump_coverage.py`, etc.).

## Datasets / inputs

The TC2214 typhoon simulation scenario over the western Pacific
(September 2022) is used as the evaluation case. It consists of
GPV-derived initial and boundary fields on an $899\times899\times128$
grid (about $10^8$ grid points) with ~2 km horizontal resolution. All
major physical processes (cloud microphysics, radiation, turbulence,
surface processes) are active. The full simulation runs for 360 time
steps ($\Delta t_{\text{big}} = 5$ s, 30 min of simulated time,
`etime=1800.e0` in `test_real/user.conf`).

The GPV input files (two binary files, ~108 MB each) are not included
in this repository due to size; they are distributed separately on
request to the AD/AE Reproducibility committee. The GPV-to-CReSS
preprocessing pipeline (`FLOW-pripro.sh`, `gridata.exe`) is part of
this repository and can regenerate the preprocessed inputs from raw
GPV data.

## Installation and deployment

1. `git clone https://github.com/hoshino-UTokyo/CReSS-GPU-SC26-AD.git`
   on a Lustre/NFS filesystem visible from the GPU node. Set
   `TOP=$(pwd)/CReSS-GPU-SC26-AD`.
2. `module load nvidia/25.9` (or equivalent) to make `nvfortran`,
   `nsys`, `ncu`, and HPC-X visible.
3. `cd $TOP && sh ispack_download.sh && ./configure.csh` builds
   ISPACK locally.
4. Compile flags are centralised in `compile.conf`; the GPU build
   uses

   ```
   -fast -mp -Mpreprocess -Mbyteswapio
   -mcmodel=medium -Minfo=accel -acc
   -gpu=managed
   ```

   as documented in Sec. IV-A. The CPU baseline build is obtained by
   removing `-acc -gpu=managed`.
5. Build the application: `cd test_real && make` (or run
   `MIYABI-solver.sh`, which invokes the build before launching the
   job).

## Workflow

The reproduction consists of three tasks:

- **$T_1$ (input staging)** Obtain GPV files (provided in the
  archive) and confirm `user.conf` parameters (`xdim=899`, `ydim=899`,
  `zdim=128`, `exprim='TC2214'`, `etime=1800.e0`).
- **$T_2$ (simulation)** Run the GPU build of CReSS for 360 steps
  using `test_real/MIYABI-solver.sh`. The MPI configuration is fixed
  to a single rank (`NODES=1, PROCS=1, CORES=72`) so that the run
  fits in a single GH200 node and matches the configuration used in
  the paper.
- **$T_3$ (validation and profiling)** Post-process
  `result/TC2214.mon.check.txt` with the validation procedure of
  Sec. III-A; profile with `nsys profile -t cuda,openacc -o run …`
  and rank kernels by accumulated GPU time. Roofline ratios are
  computed from the same `nsys` output.

Dependencies: $T_1 \rightarrow T_2 \rightarrow T_3$. No randomization
is involved; the workflow is deterministic and a single run is
sufficient. Performance numbers in Sec. IV-D are reported as the
median per-step time over the 360 steps.

## Outputs

The simulation produces `TC2214.mon.check.txt`,
`TC2214.dmp.check.txt`, and `TC2214.geography.check.txt` under
`test_real/result/`. The validation script verifies that
$|a_1|, |a_2| < 10^{-4}$. The `nsys` report yields the top-10 kernel
ranking that reproduces Fig. 1 and the bandwidth distribution used to
argue that the implementation is broadly reasonable as an OpenACC
realization. The CPU baseline run on the same node, obtained by
toggling `-acc -gpu=managed` off in `compile.conf`, provides the
9.51 s median per-step time that anchors the 5.1× GPU speedup
statement.

## Excluded content

Per-replicate Claude Code session transcripts (`.claude/` directories)
and multi-GB intermediate dump files are intentionally not included in
this repository due to size and personal-data considerations; they are
available on request to the AD/AE Reproducibility committee. The same
applies to the raw GPV input files mentioned above.

## License

See `LICENSE`. This repository is released under a paper-verification-only
license; redistribution and derivative use require explicit permission.
