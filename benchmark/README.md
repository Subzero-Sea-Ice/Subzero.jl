# Benchmarks: `timestep_sim!` on the GPU vs the CPU baseline

`timestep_sim!` is being ported to KernelAbstractions kernels one step at a time (see
`AGENTS.md`). These benchmarks measure every step of a timestep with:

- **baseline**: Subzero before the GPU work (`main` at ae96ff2), on the CPU with Julia threads.
  It doesn't change, so it runs once per machine.
- **current**: the checked-out Subzero, on the CPU backend and on NVIDIA (CUDA.jl) or AMD
  (AMDGPU.jl) GPUs. Steps that aren't ported run the same CPU code as the baseline.

Each step is a file in [`steps/`](steps/) with notes in `steps/NN_<name>.md`; results are in
`results/<machine>/`.

## Overview

The reference machine, `laptop-i7-12700H` (see [Machines](#machines)). Every machine has
its own table in `results/<machine>/README.md`.

<!-- overview:start -->
| step | ported | baseline 1 thread | baseline 12 threads | current CPU 1 thread | current CPU 12 threads | current CUDA Float64 | current CUDA Float32 | speedup | share of timestep | check | measured |
|---|:-:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|---|
| [`timestep_sim!` (whole timestep)](steps/00_timestep_sim.md) | partly | 42,990 | 24,750 | 49,860 | 21,010 | 22,520 | – | 1.18× | 100% | ✓ | 2026-10-15 `9d4fe06` |
| [`add_ghosts!` and removing the ghosts](steps/10_ghosts.md) | no | 12.1 | 15.1 | 11.6 | 16.5 | – | – | – | <1% | ✓ | 2026-10-15 `9d4fe06` |
| [`timestep_collisions!`](steps/20_collisions.md) | no | 37,130 | 8,716 | 43,080 | 8,672 | – | – | – | 35% | ✓ | 2026-10-15 `9d4fe06` |
| [`timestep_ridging_rafting!`](steps/30_ridging_rafting.md) | no | 4,054 | 4,780 | 3,925 | 5,061 | – | – | – | 19% | ✓ | 2026-10-15 `9d4fe06` |
| [`timestep_coupling!`](steps/40_coupling.md) | no | 7,736 | 8,255 | 7,996 | 8,374 | – | – | – | 33% | ✓ | 2026-10-15 `9d4fe06` |
| [`timestep_floe_properties!`](steps/50_floe_properties.md) | yes | 263 | 48.5 | 39.3 | 9.34 | 6.57 | 1.64 | 7.39× | <1% | ✓ | 2026-10-15 `9d4fe06` |
| [fracture check (`determine_fractures`, part of `fracture_floes!`)](steps/60_fracture_check.md) | no | 126 | 138 | 129 | 140 | – | – | – | 1% | ✓ | 2026-10-15 `9d4fe06` |
| [`fracture_floes!` (including the fracture check)](steps/70_fracture_floes.md) | no | 191 | 175 | 187 | 176 | – | – | – | 1% | – | 2026-10-15 `9d4fe06` |
| [`timestep_welding!`](steps/80_welding.md) | no | 791 | 702 | 821 | 680 | – | – | – | 3% | ✓ | 2026-10-15 `9d4fe06` |
| [`simplify_floes!`](steps/90_simplify.md) | no | 2.65 | 2.60 | 2.54 | 2.76 | – | – | – | <1% | ✓ | 2026-10-15 `9d4fe06` |
| [host ↔ device: `FixedWidthFloes`, `adapt`, `update_floes!`](steps/95_transfer.md) | yes | – | – | 30.8 | 46.6 | 78.0 | – | – | – | – | 2026-10-15 `9d4fe06` |

Median time of one call in ms for 100,000 floes. ≈: extrapolated from fewer floes as t ∝ nᵇ, with the exponent b.

Machine `laptop-i7-12700H`: 12th Gen Intel(R) Core(TM) i7-12700H, 1 socket, hyperthreading on, run on CPUs 0-11; GPU: NVIDIA GeForce RTX 3050 Ti Laptop GPU; the reference input `2f2e3e6e9230`.
<!-- overview:end -->

- **baseline / current**: median time of one call in ms. The GPU columns show the ported
  kernels only; the current CPU columns show the kernels on the CPU backend for ported steps,
  and the baseline code otherwise.
- **speedup**: baseline with all threads / the fastest current Float64 (ported steps only).
- **share of timestep**: baseline with all threads / the baseline of
  [`00_timestep_sim`](steps/00_timestep_sim.md), which runs every step once except welding,
  so the other shares add up to roughly 100% (ridging/rafting is overestimated, see its
  notes; `70_fracture_floes` includes `60_fracture_check`).
- **check**: the current code gives the baseline's result on every backend (✓), see
  [What is measured](#what-is-measured).
- **measured**: date and Subzero commit of the current results.

## Input

`make_input.jl` builds the docs example `docs/literate/examples/shear_flow.jl` on a larger
domain, with the same settings: periodic boundaries, shear ocean flow (0 → 0.5 → 0 m/s),
`Δgrid = 2e3`, `hmean = 0.25`, `Δh = 0`, `Δt = 20`, concentration 0.75, `Float64`,
`SubGridPointsGenerator(npoint_per_cell = 2)`, `rng = Xoshiro(1)`; collisions, coupling and
simplification on, fractures, ridging/rafting and welding off. Instead of 50 floes on
1e5 × 1e5 m, it asks for `nfloes = 10_000` on 1.416e6 × 1.416e6 m (same density, so the
same floe size) and gets **10,029 floes** with 7.0 vertices on average. It then runs 100
timesteps so that velocities, forces and interactions are realistic (2.4 interactions per floe
on average, at most 7), and saves the floe fields as plain arrays in `floes_10k.jld2`, so
both Subzero versions load exactly the same floes.

**Reproducibility.** `run.sh` makes the input with 1 thread and the baseline, and
`make_input.jl` seeds the global rng, so a machine makes the same file every time. With more
threads, the spin-up depends on the order in which threads get a lock in
`timestep_collisions!`; with the current code, the input would change with every commit that
changes round-off; and before Subzero passes its `rng` to `GO.voronoi`, the tesselation uses
the global rng. Other package versions or CPUs can still change the last bits, which the
spin-up amplifies. So every result records the first 12 hex digits of the input's SHA-256
(`input`), and each machine's table says if it isn't the reference input (`2f2e3e6e9230`).
The file is 90 MB and not committed: make it with `run.sh` (or ask for it), check
`sha256sum floes_10k.jld2`, and copy it to the work directory of other machines rather than
making it there.

**Size.** Typical simulations have 100,000 to 1,000,000 floes. The default is **100,000
floes**: the 10,029-floe field copied 10 times, each shifted by `Lx` in x, with the grid,
periodic domain and shear flow copied the same way (`make_model` in `lib/input.jl`), so the
density stays realistic. Set `BENCH_SIZES` for other sizes, e.g. 1,000,000 on an HPC system
(see [Memory use](#memory-use)).

## What is measured

BenchmarkTools, 20 samples per benchmark (`evals = 1`, at most 60 s). The garbage collector
runs before every sample, and setup (e.g. a `deepcopy` of floes that a step changes) isn't
timed. Every step file says what it measures and what its setup is. The baseline's per-event
`@info` logging is silenced with a `NullLogger`.

**Extrapolation.** Collisions check every pair of floes (O(n²)), so some baseline steps are
too slow for 1,000,000 floes. They set `extrapolate_above` (100,000) and `extrapolate_from`
(25,000, 50,000, 100,000): above that size, they run on the smaller sizes, fit `t = a nᵇ` to
the median times, and the overview shows `≈` with `b` in a footnote. The same applies to the
current code of unported steps. Ported code must run at the sizes in `BENCH_SIZES`.

**Correctness check.** The baseline saves its result on 1,000 floes (or `check_n`) in the
work directory. The current code compares its result on every backend with it: the maximum
difference relative to the largest value of each field, and whether all fields are
`isapprox`. Steps that use the rng in a threaded loop can't be checked.

## Machines

Results of different machines can't be compared, so each has its own folder
`results/<machine>/` (`BENCH_MACHINE`, default: the host name) with its overview table in
`README.md`. The results files record the date, Subzero commit, CPU, the CPUs used
(`BENCH_CPUS` or the Slurm job's), hyperthreading, the number of sockets of those CPUs and
the GPU.

The reference machine (`REFERENCE` in `summarize.jl`), `laptop-i7-12700H`:

- Laptop, on AC power, `balanced` platform profile, `powersave` CPU governor.
- Intel Core i7-12700H (6 P-cores with hyperthreading + 8 E-cores, 20 threads), 32 GB RAM.
  The benchmarks use the 12 threads of the P-cores only (CPUs 0–11), see
  [CPU cores](#cpu-cores).
- NVIDIA GeForce RTX 3050 Ti Laptop GPU, 4 GB, driver 580.178.04. Its Float64 throughput is
  1/64 of Float32, as usual for consumer GPUs, so check correctness in Float64 and look at
  Float32 for speed. A datacenter GPU (A100/H100) with full-rate Float64 performs very
  differently.
- Julia 1.12.7, CUDA.jl 6.4.2 (CUDA runtime 13.4), KernelAbstractions 0.9.45,
  BenchmarkTools 1.8.0.

### CPU cores

The i7-12700H has P-cores (up to 4.7 GHz, hyperthreading, CPUs 0–11) and E-cores (up to
3.5 GHz, no hyperthreading, CPUs 12–19). Julia doesn't pin threads, so an unpinned 1-thread
run lands on either type, and a run on all cores mixes them. Median ms for 50,000 floes
(2026-10-15, pinned with `taskset`):

| CPUs | collisions (baseline) | coupling (baseline) | floe properties (baseline) | floe properties (current CPU kernels) |
|---|--:|--:|--:|--:|
| 1 thread, unpinned | 7,669 | 2,640 | 131 | 19.5 |
| 1 thread, P-core | 8,466 | 2,500 | 121 | 17.5 |
| 1 thread, E-core | 11,963 | 3,890 | 200 | 33.4 |
| 6 threads, P-cores, no hyperthreading | 2,710 | 2,686 | 26.7 | 6.73 |
| 12 threads, P-cores with hyperthreading | 2,087 | 2,795 | 24.3 | 3.63 |
| 14 threads, P-cores without hyperthreading + E-cores | 1,766 | 2,760 | 20.5 | 3.59 |
| 20 threads, unpinned | 1,611 | 2,640 | 16.8 | 2.60 |

An E-core is 1.4–1.9× slower than a P-core; hyperthreading makes 12 threads 1.3–1.9× faster
than 6; the E-cores add another 10–40% but mix two speeds; coupling doesn't use threads. So
`run.sh` runs on the P-cores only (`taskset -c 0-11`, 12 threads), also for 1 thread: that's
comparable between runs and closer to an HPC node with one type of core. Set `BENCH_CPUS` to
choose, e.g. `BENCH_CPUS=0-19` for all cores.

On other machines:

- **One type of core** (most desktops and HPC nodes): all CPUs the process may use (in Slurm,
  the job's), one thread each, including hyperthreads. The warnings below only look at those
  CPUs.
- **Other fast + slow core CPUs** (AMD Zen 5 + Zen 5c, ARM big.LITTLE) aren't detected:
  `run.sh` warns if maximum frequencies differ by more than 20%; set `BENCH_CPUS` to the fast
  ones.
- **More than one socket**: unpinned threads, and their memory, move between sockets.
  `run.sh` warns; set `BENCH_CPUS` to one socket's CPUs (see `lscpu`) for less noise.
- **macOS**: no `taskset`, so no pinning and no `BENCH_CPUS`; hyperthreading and sockets are
  recorded as unknown.

## How to run

```sh
REPO=/path/to/Subzero.jl            # checked out at the commit to benchmark
WORK=$(mktemp -d)                   # work directory for environments, input, and logs
$REPO/benchmark/setup.sh "$WORK"    # worktree of the baseline + env-baseline and env-current
cd "$WORK" && BENCH_MACHINE=<name> $REPO/benchmark/run.sh   # e.g. laptop-i7-12700H
git -C $REPO worktree remove "$WORK/subzero-baseline"   # clean up afterwards
```

On Snellius, `sbatch benchmark/snellius.sh` does all of this on a quarter of a `gpu_a100`
node (18 cores, one A100); see the comments in the script.

`setup.sh` adds CUDA.jl to the current environment if `nvidia-smi` is found and AMDGPU.jl if
`rocminfo` or `rocm-smi` is; choose with `BENCH_GPU` (e.g. `BENCH_GPU=AMDGPU`, or `""` for
none). Its second argument is the baseline commit (default ae96ff2).

`run.sh`, on the CPUs in `BENCH_CPUS` with `BENCH_THREADS` threads:

1. `make_input.jl 10000 100 floes_10k.jld2`: the input, with 1 thread and the baseline (a few
   minutes), unless `floes_10k.jld2` is in the work directory.
2. `bench.jl baseline` with 1 and `BENCH_THREADS` threads, skipping steps with baseline
   results for that number of threads and this input checksum.
3. `bench.jl current` with 1 thread (`BENCH_BACKENDS=CPU`) and `BENCH_THREADS` threads (CPU
   and GPUs).
4. `summarize.jl` updates the overview tables.

Results go to `results/<machine>/NN_<name>.toml` here, output to `logs/` in the work
directory. To run some steps only, e.g. after porting one: `BENCH_STEPS=collisions run.sh`.

| Variable | Meaning | Default |
|---|---|---|
| `BENCH_MACHINE` | machine name, for the results folder | the host name |
| `BENCH_SIZES` | numbers of floes | `100000` |
| `BENCH_STEPS` | steps to run, with or without their number | all |
| `BENCH_SAMPLES` | samples per benchmark | 20 |
| `BENCH_SECONDS` | maximum time per benchmark in s | 60 |
| `BENCH_BACKENDS` | from `CPU`, `CUDA`, `AMDGPU` | `CPU` + each of CUDA.jl and AMDGPU.jl in the environment with a working GPU |
| `BENCH_FLOES` | input file | `floes_10k.jld2` |
| `BENCH_FORCE=1` | rerun the baseline even if its results exist | |
| `BENCH_CPUS` (`run.sh`) | CPUs to run on, as a `taskset` list | the fast cores if there are fast and slow ones, else all |
| `BENCH_THREADS` (`run.sh`) | number of threads | the number of CPUs in `BENCH_CPUS` |

Lists are comma-separated.

## Adding or porting a step

A step is a file `steps/NN_<name>.jl` that calls `register_step!` (see `Step` in
`lib/harness.jl`). `bench.jl` runs the files in name order, so there's no list to edit. The
numbers follow `timestep_sim!`, with gaps so that a step can be split later.

```julia
register_step!(;
    title = "`timestep_floe_properties!`",   # for the overview table
    baseline = function (ctx)                # baseline code, CPU
        floes = load_floes(ctx.n)
        trial = @bench(Subzero.timestep_floe_properties!(f, 1, 20, $fs), f = deepcopy($floes))
        record!(ctx, "timestep_floe_properties!", trial; primary = true)
    end,
    current = function (ctx, backend, FT)    # current code; leave out until the step is ported
        ...
        record!(ctx, "kernels", trial; FT, primary = true)
    end,
    baseline_result = n -> Dict(...),        # arrays after one call, for the correctness check
    current_result = (n, backend) -> Dict(...),
    float32 = true,                          # also run `current` in Float32 on GPUs
    extrapolate_above = typemax(Int),        # if the baseline is slow: e.g. 100_000, and
    extrapolate_from = Int[],                # e.g. [25_000, 50_000, 100_000]
)
```

- `record!` keeps a measurement; the one with `primary = true` goes in the overview table,
  others (e.g. parts of the step) only in the results file.
- Step files are loaded in both environments, so functions that only exist in the current
  Subzero can be used in function bodies, not at the top level.
- `lib/input.jl` has the input (`load_floes(n)`, `make_model(floes, n)`, `ghosted_floes(n)`,
  `collided_floes(n)`, `input_simulation(n; kwargs...)`); `lib/backends.jl` has helpers for
  floes on the device (`to_float`, `reset!`, `free!`).

To port a step, add `current` and `current_result`, run `BENCH_STEPS=<step> run.sh`, and write
what you see in `steps/NN_<name>.md`.

### Keeping merge conflicts small

Several steps are ported at the same time, so:

1. A porting PR changes only `steps/NN_<step>.{jl,md}`, `results/<machine>/NN_<step>.toml`
   and that step's row in the overview tables. Helpers start in the step file and move to
   `lib/` in a separate small PR if other steps need them.
2. The tables have a row for every step from the start, so PRs only change rows. Neighbouring
   rows still conflict in git, but the results files never do, so take either side of the
   tables and make them again (`summarize.jl --check` fails if they disagree):

   ```sh
   git checkout --theirs benchmark/README.md benchmark/results/*/README.md
   julia benchmark/summarize.jl && git add benchmark/README.md benchmark/results/*/README.md
   ```

3. The baseline is pinned to ae96ff2, not `main`, which gets the GPU code. Its results are
   committed, so a PR only runs the current code of its step on a machine with baseline
   results.
4. Without the reference machine, run `BENCH_MACHINE=<name> BENCH_STEPS=<step> run.sh`, which
   also runs that step's baseline there (5–30 minutes), and commit `results/<name>/`. The
   reference row is updated later, e.g. in a refresh of all steps.
5. Dependencies of the benchmarks go in `Project.toml` in this folder.

## Memory use

Measured with an earlier version (`timestep_floe_properties!` only, 20 threads, CPU and
CUDA backends):

| Floes | host (max RSS) | GPU |
|---:|---:|---:|
| 100,000 | 2.7 GB | 0.32 GB |
| 200,000 | 3.8 GB | 0.49 GB |
| 400,000 | 5.7 GB | 0.86 GB |

That's about 1 GB (host) and 0.18 GB (GPU) per 100,000 floes. Steps that need the grid, ocean
and atmosphere need more: for 1,000,000 floes, the grid has about 50 million cells, so
run 1,000,000 floes on an HPC system. `bench.jl` prints its peak memory use at the end.
