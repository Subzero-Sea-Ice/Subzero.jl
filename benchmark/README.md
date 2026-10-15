# Benchmarks: `timestep_sim!` on the GPU vs the CPU baseline

`timestep_sim!` is being ported to KernelAbstractions kernels one step at a time (see
`AGENTS.md`). These benchmarks measure every step of a timestep twice:

- **baseline**: Subzero before the GPU work (`main` at ae96ff2), on the CPU with Julia threads.
  It doesn't change, so it only needs to run once per machine.
- **current**: the checked-out Subzero, on the CPU backend and on NVIDIA (CUDA.jl) or AMD
  (AMDGPU.jl) GPUs. Steps that aren't ported yet run the same CPU code as the baseline.

Each step is a file in [`steps/`](steps/), with its results in `results/<machine>/` and
notes and observations in `steps/NN_<name>.md`.

## Overview

The results of the reference machine, `laptop-i7-12700H` (see
[Machines](#machines)). Every machine has its own table in `results/<machine>/README.md`.

<!-- overview:start -->
| step | ported | baseline 1 thread | baseline 12 threads | current CPU 1 thread | current CPU 12 threads | current CUDA Float64 | current CUDA Float32 | speedup | share of timestep | check | measured |
|---|:-:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|---|
| [`timestep_sim!` (whole timestep)](steps/00_timestep_sim.md) | partly | 45,260 | 17,970 | 41,910 | 21,880 | 18,120 | – | 0.992× | 100% | ✓ | 2026-10-15 `1160b6a` |
| [`add_ghosts!` and removing the ghosts](steps/10_ghosts.md) | no | 12.6 | 16.0 | 11.9 | 16.3 | – | – | – | <1% | ✓ | 2026-10-15 `1160b6a` |
| [`timestep_collisions!`](steps/20_collisions.md) | no | 35,900 | 9,551 | 34,480 | 12,030 | – | – | – | 53% | ✓ | 2026-10-15 `1160b6a` |
| [`timestep_ridging_rafting!`](steps/30_ridging_rafting.md) | no | 2,845 | 4,554 | 2,726 | 4,571 | – | – | – | 25% | ✓ | 2026-10-15 `1160b6a` |
| [`timestep_coupling!`](steps/40_coupling.md) | no | 6,156 | 5,399 | 5,488 | 5,211 | – | – | – | 30% | ✓ | 2026-10-15 `1160b6a` |
| [`timestep_floe_properties!`](steps/50_floe_properties.md) | yes | 383 | 49.0 | 40.8 | 6.20 | 6.58 | 1.64 | 7.91× | <1% | ✓ | 2026-10-15 `1160b6a` |
| [fracture check (`determine_fractures`, part of `fracture_floes!`)](steps/60_fracture_check.md) | no | 127 | 142 | 132 | 154 | – | – | – | 1% | ✓ | 2026-10-15 `1160b6a` |
| [`fracture_floes!` (including the fracture check)](steps/70_fracture_floes.md) | no | 336 | 178 | 202 | 212 | – | – | – | 1% | – | 2026-10-15 `1160b6a` |
| [`timestep_welding!`](steps/80_welding.md) | no | 662 | 705 | 660 | 701 | – | – | – | 4% | ✓ | 2026-10-15 `1160b6a` |
| [`simplify_floes!`](steps/90_simplify.md) | no | 2.48 | 2.60 | 2.48 | 2.46 | – | – | – | <1% | ✓ | 2026-10-15 `1160b6a` |
| [host ↔ device: `FixedWidthFloes`, `adapt`, `update_floes!`](steps/95_transfer.md) | yes | – | – | 31.3 | 34.5 | 63.1 | – | – | – | – | 2026-10-15 `1160b6a` |

Median time of one call in ms for 100,000 floes. ≈: extrapolated from fewer floes as t ∝ nᵇ, with the exponent b.

Machine `laptop-i7-12700H`: 12th Gen Intel(R) Core(TM) i7-12700H, 1 socket, hyperthreading on, run on CPUs 0-11; GPU: NVIDIA GeForce RTX 3050 Ti Laptop GPU; the reference input `d1805bc093a4`.
<!-- overview:end -->

- **baseline / current**: median time of one call in ms. The current Float64 and Float32
  columns for GPUs show the ported kernels only; the current CPU columns show the kernels on
  the CPU backend for ported steps, and the baseline code for steps that aren't ported.
- **speedup**: baseline with all threads / the fastest current Float64 (ported steps only).
- **share of timestep**: baseline with all threads / the baseline of
  [`00_timestep_sim`](steps/00_timestep_sim.md). That runs every step once, except welding, so
  the shares of the other steps add up to roughly 100% (ridging/rafting is overestimated, see
  its notes). `70_fracture_floes` includes
  `60_fracture_check`.
- **check**: the current code gives the same result as the baseline on every backend (✓), see
  [What is measured](#what-is-measured).
- **measured**: date and Subzero commit of the current results.

## Input

`make_input.jl` builds the docs example `docs/literate/examples/shear_flow.jl` on a larger
domain:

- Same settings as the example: periodic boundaries, shear ocean flow (0 → 0.5 → 0 m/s),
  `Δgrid = 2e3`, `hmean = 0.25`, `Δh = 0`, `Δt = 20`, concentration 0.75, `Float64`,
  `SubGridPointsGenerator(npoint_per_cell = 2)`, `rng = Xoshiro(1)`. Collisions, coupling
  and simplification are on; fractures, ridging/rafting and welding are off.
- Instead of 50 floes on a 1e5 × 1e5 m domain, `nfloes = 10_000` on a 1.416e6 × 1.416e6 m
  domain (same floe density, so the floes have the same size). The generator made
  **10,029 floes**, with 7.0 vertices per floe on average.
- The simulation runs 100 timesteps so that velocities, forces, and collision interactions
  are realistic. After that, floes have 2.4 interactions on average (at most 7).
- The floe fields are saved as plain arrays in `floes_10k.jld2`, so both Subzero versions
  load exactly the same floes.

**Making the same input again.** `run.sh` runs `make_input.jl` with 1 thread and the
baseline, and `make_input.jl` seeds the global rng. Then the same machine makes the same file
every time:

- With more threads, the spin-up depends on the order in which the threads get a lock in
  `timestep_collisions!`, and the 100 timesteps amplify that.
- With the current code, the input would change with every commit that changes round-off.
- Before Subzero passes its `rng` to `GO.voronoi`, the Voronoi tesselation uses the global
  rng.

Other package versions or CPUs can still change the last bits, and the spin-up amplifies
them. So every result records the first 12 digits of the SHA-256 of the input file
(`input`), and each machine's table says if it isn't the input of the reference results
(`d1805bc093a4`). That input was made before these changes (with 20 threads and the current
code), so it can't be made again: it's 24 MB and not committed, so ask for it, and copy it
to the work directory of every machine.

Typical simulations have 100,000 to 1,000,000 floes. The benchmarks use **100,000 floes** by
default, made by copying the 10,029-floe field 10 times, each copy shifted by `Lx` in x. The
grid, the periodic domain and the ocean's shear flow are copied in the same way
(`make_model` in `lib/input.jl`), so the floe density stays realistic. Set `BENCH_SIZES` for
other sizes, e.g. 1,000,000 floes on an HPC system (see [Memory use](#memory-use)).

## What is measured

All numbers come from BenchmarkTools, with 20 samples per benchmark (`evals = 1`, at most
60 s per benchmark). The garbage collector runs before every sample, which isn't timed.
Setup, such as making a `deepcopy` of floes that a step changes, isn't timed either. Every
step file says what it measures and what its setup is. The baseline's per-event `@info`
logging is silenced with a `NullLogger`.

**Too slow for more floes.** Some baseline steps get too slow for 1,000,000 floes: collisions
check every pair of floes, so their time grows with the square of the number of floes. Those
steps set `extrapolate_above` (100,000) and `extrapolate_from` (25,000, 50,000 and 100,000
floes). For more floes than `extrapolate_above`, the baseline runs on the `extrapolate_from`
sizes, `t = a nᵇ` is fitted to their median times, and the overview shows the extrapolated
time as `≈`, with the exponent `b` in a footnote. The same applies to the current code of a
step that isn't ported yet. Ported code must run at the sizes in `BENCH_SIZES`.

**Correctness check.** The baseline saves its result on 1,000 floes (or `check_n`) in the
work directory, and the current code compares its result on every backend with it: the
maximum difference relative to the largest value of each field, and whether all fields are
`isapprox`. Steps that use the rng in a threaded loop can't be checked.

## Machines

Results of different machines can't be compared, so every machine has its own folder,
`results/<machine>/`, named by `BENCH_MACHINE` (default: the host name). Its results files
record the date, Subzero commit, CPU, the CPUs the benchmark ran on (set by `BENCH_CPUS`, or
those of a Slurm job), whether hyperthreading was on, the number of sockets of those CPUs and
the GPU, and its `README.md` has its overview table.

The table above is the reference machine, `laptop-i7-12700H` (`REFERENCE` in
`summarize.jl`):

- Laptop, on AC power, `balanced` platform profile, `powersave` CPU governor
- CPU: Intel Core i7-12700H (6 P-cores with hyperthreading + 8 E-cores, 20 threads), 32 GB
  RAM. The benchmarks run on the 12 threads of the P-cores only (CPUs 0–11), see
  [CPU cores](#cpu-cores).
- GPU: NVIDIA GeForce RTX 3050 Ti Laptop GPU, 4 GB, driver 580.178.04. Consumer GPUs have
  low Float64 throughput (1/64 of Float32 for this one), so check correctness in Float64 and
  look at Float32 for speed. A datacenter GPU (A100/H100) with full-rate Float64 would
  perform very differently.
- Julia 1.12.7, CUDA.jl 6.4.2 (CUDA runtime 13.4), KernelAbstractions 0.9.45,
  BenchmarkTools 1.8.0

### CPU cores

This CPU has fast cores (P-cores, up to 4.7 GHz, with hyperthreading: CPUs 0–11) and slow
cores (E-cores, up to 3.5 GHz, no hyperthreading: CPUs 12–19). Julia doesn't pin its threads
to cores, so without pinning, a 1-thread run can land on either type, and a run on all
cores mixes them. Median times in ms for 50,000 floes (2026-10-15, `taskset` for pinning):

| CPUs | collisions (baseline) | coupling (baseline) | floe properties (baseline) | floe properties (current CPU kernels) |
|---|--:|--:|--:|--:|
| 1 thread, unpinned | 7,669 | 2,640 | 131 | 19.5 |
| 1 thread, P-core | 8,466 | 2,500 | 121 | 17.5 |
| 1 thread, E-core | 11,963 | 3,890 | 200 | 33.4 |
| 6 threads, P-cores, no hyperthreading | 2,710 | 2,686 | 26.7 | 6.73 |
| 12 threads, P-cores with hyperthreading | 2,087 | 2,795 | 24.3 | 3.63 |
| 14 threads, P-cores without hyperthreading + E-cores | 1,766 | 2,760 | 20.5 | 3.59 |
| 20 threads, unpinned | 1,611 | 2,640 | 16.8 | 2.60 |

- An E-core is 1.4–1.9× slower than a P-core.
- Hyperthreading helps: 12 threads on the 6 P-cores are 1.3–1.9× faster than 6.
- The E-cores add another 10–40%, but mix two speeds of cores.
- Coupling doesn't use threads.

So `run.sh` runs Julia on the P-cores only (`taskset -c 0-11`, 12 threads), also for the
1-thread runs: comparable from run to run, and closer to an HPC node, which has one type of
core. Set `BENCH_CPUS` to choose, e.g. `BENCH_CPUS=0-19` for all cores of this CPU.

On other machines:

- **One type of core** (most desktops and HPC nodes): `run.sh` uses all CPUs that the process
  may use, with one thread per CPU, including hyperthreads if they're on. In a Slurm job,
  that's the CPUs of the job. The warnings below only look at those CPUs.
- **Other CPUs with fast and slow cores** (e.g. AMD Zen 5 + Zen 5c, ARM big.LITTLE) aren't
  detected. `run.sh` warns if the maximum frequencies of the CPUs differ by more than 20%;
  then set `BENCH_CPUS` to the fast ones.
- **CPUs on more than one socket**: threads that aren't pinned move between the sockets, and their
  memory with them. `run.sh` warns; set `BENCH_CPUS` to the CPUs of one socket (see `lscpu`)
  for less noise.
- **macOS**: there's no `taskset`, so `run.sh` doesn't pin and `BENCH_CPUS` doesn't work.
  Hyperthreading and the number of sockets are recorded as unknown.

## How to run

```sh
REPO=/path/to/Subzero.jl            # checked out at the commit to benchmark
WORK=$(mktemp -d)                   # work directory for environments, input, and logs
$REPO/benchmark/setup.sh "$WORK"    # worktree of the baseline + env-baseline and env-current
cd "$WORK" && BENCH_MACHINE=<name> $REPO/benchmark/run.sh   # e.g. laptop-i7-12700H
git -C $REPO worktree remove "$WORK/subzero-baseline"   # clean up afterwards
```

On Snellius, `sbatch benchmark/snellius.sh` runs all of this in a Slurm job, on a quarter of a
`gpu_a100` node (18 cores and one A100); see the comments in the script.

`run.sh` does the following, on the CPUs in `BENCH_CPUS` with `BENCH_THREADS` threads (see
[CPU cores](#cpu-cores)):

1. `make_input.jl 10000 100 floes_10k.jld2` generates the input with 1 thread and the
   baseline (a few minutes). This is skipped if `floes_10k.jld2` exists in the work
   directory.
2. `bench.jl baseline` with 1 and `BENCH_THREADS` threads. Steps whose baseline results exist
   for that number of threads are skipped.
3. `bench.jl current` with 1 thread (CPU only, `BENCH_BACKENDS=CPU`) and `BENCH_THREADS`
   threads (CPU and GPUs).
4. `summarize.jl` updates the overview tables.

The results go to `results/<machine>/NN_<name>.toml` in this folder, the output to `logs/` in
the work directory. To run only some steps, e.g. after porting one:
`BENCH_STEPS=collisions run.sh` (names with or without their number, comma-separated).

`setup.sh` adds CUDA.jl to the current environment if `nvidia-smi` is found, and AMDGPU.jl
if `rocminfo` or `rocm-smi` is found. To choose yourself, set `BENCH_GPU`, for example
`BENCH_GPU=AMDGPU setup.sh "$WORK"` (or `BENCH_GPU=""` for no GPU package). The second
argument of `setup.sh` sets the baseline commit (default ae96ff2).

Environment variables for `run.sh` and `bench.jl`:

- `BENCH_MACHINE`: name of the machine, for its results folder (default: the host name)
- `BENCH_SIZES`: comma-separated numbers of floes (default `100000`)
- `BENCH_STEPS`: comma-separated steps to run (default: all)
- `BENCH_SAMPLES`: samples per benchmark (default 20)
- `BENCH_SECONDS`: maximum time per benchmark in s (default 60)
- `BENCH_BACKENDS`: comma-separated, from `CPU`, `CUDA`, `AMDGPU` (default `CPU` plus each of
  CUDA.jl and AMDGPU.jl that is in the environment and has a working GPU)
- `BENCH_FLOES`: input file (default `floes_10k.jld2`)
- `BENCH_FORCE=1`: rerun the baseline even if its results exist
- `BENCH_CPUS` (`run.sh` only): CPUs to run on, as a `taskset` list (default: the fast cores
  of a CPU with fast and slow cores, otherwise all CPUs)
- `BENCH_THREADS` (`run.sh` only): number of threads (default: the number of CPUs in
  `BENCH_CPUS`)

## Adding or porting a step

A step is a file `steps/NN_<name>.jl` that calls `register_step!` (see `Step` in
`lib/harness.jl`). `bench.jl` finds the files itself and runs them in the order of their
names, so there's no list of steps to edit. The numbers follow the order of `timestep_sim!`,
with gaps so that a step can be split later.

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

- `record!` keeps a measurement. The one with `primary = true` goes in the overview table;
  others (e.g. parts of the step) are only in the results file.
- The step files are loaded in both environments. Functions that only exist in the current
  Subzero are fine inside function bodies, but not at the top level of a step file.
- `lib/input.jl` has the input: `load_floes(n)`, `make_model(floes, n)`, `ghosted_floes(n)`,
  `collided_floes(n)`, `input_simulation(n; kwargs...)`. `lib/backends.jl` has helpers for
  floes on the device: `to_float`, `reset!`, `free!`.

To port a step, give it a `current` function (and `current_result`), run
`BENCH_STEPS=<step> run.sh`, and write what you see in `steps/NN_<name>.md`.

### Keeping merge conflicts small

Several steps are ported at the same time, so:

1. A PR that ports a step changes only `steps/NN_<step>.{jl,md}`,
   `results/<machine>/NN_<step>.toml` and that step's row in the overview tables (made by
   `summarize.jl`). Helpers start in the step file. If other steps need them too, they move
   to `lib/` in a separate small PR.
2. The overview tables have a row for every step from the start, so PRs only change rows.
   Two PRs that change neighbouring rows still conflict, because git treats changes to
   neighbouring lines as one change. The results files never conflict, so take either side
   of the tables and make them again:

   ```sh
   git checkout --theirs benchmark/README.md benchmark/results/*/README.md
   julia benchmark/summarize.jl && git add benchmark/README.md benchmark/results/*/README.md
   ```

   `julia benchmark/summarize.jl --check` fails if a table and the results files disagree.
3. The baseline is pinned to ae96ff2, not `main`, because `main` gets the GPU code. Its
   results are in the results files, so a PR only runs the current code of its own step on a
   machine that has baseline results.
4. Without the reference machine, a PR adds the results of its step on another machine: run
   `BENCH_MACHINE=<name> BENCH_STEPS=<step> run.sh`, which also runs the baseline of that
   step on that machine (5–30 minutes per step), and commit `results/<name>/`. The reference
   machine's row is updated later, e.g. in a refresh of all steps.
5. Dependencies of the benchmarks go in `Project.toml` in this folder.

## Memory use

Measured with the earlier version of these benchmarks (`timestep_floe_properties!` only),
with 20 threads and the CPU and CUDA backends:

| Floes | host (max RSS) | GPU |
|---:|---:|---:|
| 100,000 | 2.7 GB | 0.32 GB |
| 200,000 | 3.8 GB | 0.49 GB |
| 400,000 | 5.7 GB | 0.86 GB |

Both grow by about 1 GB (host) and 0.18 GB (GPU) per 100,000 floes. The steps that need the
grid, ocean and atmosphere need more: for 1,000,000 floes, the grid has about 50 million
cells. Run 1,000,000 floes on an HPC system. `bench.jl` prints its peak memory use at the end.
