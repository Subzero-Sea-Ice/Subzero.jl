# Benchmark: `timestep_floe_properties!` on the GPU vs `main`

Compares `timestep_floe_properties!` on the `gpu-version` branch (KernelAbstractions
kernels on `FixedWidthFloes`) with the same function on `main` (a `Threads.@threads`
loop over the floe `StructArray`). The branch is benchmarked on the CPU backend and on
NVIDIA (CUDA.jl) or AMD (AMDGPU.jl) GPUs. The results below are from an NVIDIA GPU.

## Input

`make_floes.jl` builds the docs example `docs/literate/examples/shear_flow.jl` on a larger
domain:

- Same settings as the example: periodic boundaries, shear ocean flow (0 → 0.5 → 0 m/s),
  `Δgrid = 2e3`, `hmean = 0.25`, `Δh = 0`, `Δt = 20`, concentration 0.75, `Float64`,
  `SubGridPointsGenerator(npoint_per_cell = 2)`, `rng = Xoshiro(1)`.
- Instead of 50 floes on a 1e5 × 1e5 m domain, `nfloes = 10_000` on a 1.416e6 × 1.416e6 m
  domain (same floe density, so the floes have the same size). The generator made
  **10,029 floes**, with 7.0 vertices per floe on average.
- The simulation runs 100 timesteps (on the branch, CPU backend) so that velocities,
  forces, and collision interactions are realistic. After that, floes have 2.4
  interactions on average (at most 7).
- The floe fields are saved as plain arrays in `floes_10k.jld2`, so both Subzero
  versions load exactly the same floes. A Voronoi floe field isn't reproducible with a
  fixed `rng`, so a rerun gives a slightly different field.

Sizes benchmarked: 100,000 floes by default, made by copying the 10,029-floe field 10 times,
each copy shifted by `Lx` in x. Set `BENCH_SIZES` for other sizes, e.g. 1,000,000 floes on
an HPC system (see [Memory use](#memory-use)). The function works on each floe
independently, so tiling only scales the amount of work, and keeps the floe density
realistic.

The floe settings are `FloeSettings()` (the default `DecayAreaScaledCalculator`,
`max_floe_height = 10`, `maximum_ξ = 1e-5`, as in the example).

## What is measured

All numbers come from BenchmarkTools, with 20 samples per benchmark (`evals = 1`, at most
60 s per benchmark). The garbage collector runs before every sample, which isn't timed.

- **main**: `timestep_floe_properties!(floes, 1, Δt, floe_settings)`. Every sample starts
  from a `deepcopy` of the floes, made in `setup` and not timed. Its per-event `@info`
  logging is silenced with a `NullLogger`.
- **branch, kernels only**: `timestep_floe_properties!(dev_floes, Δt, floe_settings; backend)`
  followed by `KernelAbstractions.synchronize(backend)`, with `dev_floes` already on the
  device. Every sample starts by copying the input into `dev_floes` on the device (not
  timed).
- **branch, Float32 kernels only** (GPU only): the same, with every `Float64` array of the
  `FixedWidthFloes` converted to `Float32` and `FloeSettings(Float32)`. The kernels still do
  some `Float64` maths because of `Float64` literals (`1.5Δt`), see `AGENTS.md`.
- **branch, full step**: what `run!` does around the call: `FixedWidthFloes(floes)`,
  `adapt(backend, …)`, the kernels, `synchronize`, `adapt(Array, …)`, `update_floes!`,
  and `_log_flags`. It is also measured in parts. Every sample starts by resetting the
  floes with `update_floes!` (not timed).

Correctness check: `bench_main.jl` saves the result of one call on 1,000 floes, and
`bench_branch.jl` compares the branch's result with it (centroid, height, mass, moment, α,
velocities, previous-timestep values, stresses, strain). Both CPU and CUDA match `main`,
with a maximum relative difference of 6e-15.

## Machine and versions

- Laptop, on AC power, `balanced` platform profile, `powersave` CPU governor
- CPU: Intel Core i7-12700H (6 P-cores + 8 E-cores, 20 threads), 32 GB RAM
- GPU: NVIDIA GeForce RTX 3050 Ti Laptop GPU, 4 GB, driver 580.178.04. Consumer GPUs have
  low Float64 throughput (1/64 of Float32 for this one).
- Julia 1.12.7, CUDA.jl 6.4.2 (CUDA runtime 13.4), KernelAbstractions 0.9.44,
  BenchmarkTools 1.8.0
- Subzero: `main` at ae96ff2, `gpu-version` at c717bc9 with `timestep_floe_properties!` as
  a single kernel
- Run on 2026-10-09

## Results

Median time per call in ms for 100,000 floes (Float64 unless noted). The full output is
under [Example output](#example-output).

| main, 1 thread | main, 20 threads | branch CPU, 1 thread, kernels | branch CPU, 20 threads, kernels | branch CUDA, kernels | branch CUDA Float32, kernels |
|---:|---:|---:|---:|---:|---:|
| 371 | 113 | 37.6 | 5.19 | 6.68 | 1.76 |

The full step as in `run!`, with the time spent in each part (branch, 20 threads):

| backend | full step | `FixedWidthFloes` | to device | kernels | to host | `update_floes!` |
|---|---:|---:|---:|---:|---:|---:|
| CPU  | 54.9 | 66.3 | 0    | 5.19 | 0    | 19.9 |
| CUDA | 88.8 | 52.8 | 12.3 | 6.68 | 25.5 | 20.0 |

The CPU full step is shorter than the sum of its parts. The parts are measured separately
and vary by about ±20% between runs with 20 threads.

## Observations

- **The kernels are 17–22× faster than `main`** with 20 threads, both on CUDA and on the CPU
  backend. With 1 thread, the CPU kernels are 10× faster than `main`, because they don't
  allocate.
- **Float64 is slow on this GPU.** In Float32, the CUDA kernels are 3.8× faster. A
  datacenter GPU (A100/H100) with full-rate Float64 would perform very differently.
- **In `run!`, the conversion costs much more than the kernels.** `FixedWidthFloes(floes)`
  and `update_floes!` take 73–86 ms, and copying to and from the GPU another 38 ms, against
  5–7 ms for the kernels. The GPU only pays off once more of `timestep_sim!` runs on the
  device and the floes stay there between steps.

## Memory use

Peak memory use of `bench_branch.jl` with 20 threads and the CPU and CUDA backends:

| Floes | host (max RSS) | GPU |
|---:|---:|---:|
| 100,000 | 2.7 GB | 0.32 GB |
| 200,000 | 3.8 GB | 0.49 GB |
| 400,000 | 5.7 GB | 0.86 GB |

Both grow by about 1 GB (host) and 0.18 GB (GPU) per 100,000 floes, so 1,000,000 floes need
about 12 GB of host memory and 2 GB of GPU memory. That's too much for a 32 GB laptop with a
desktop and an IDE running, so run 1,000,000 floes on an HPC system. The scripts print their
peak memory use at the end.

## Example output

The output of `run_all.sh` (the files in `results/`) for the results above.

<details>
<summary><code>main</code>, 1 thread: <code>results/main_t1.txt</code></summary>

```
main, Julia threads = 1
main: timestep_floe_properties!               n=100000  median    370.708 ms  min    356.236 ms  (3707.1 ns/floe, 14 samples)
Peak host memory (max RSS): 1.39 GB
```

</details>

<details>
<summary><code>main</code>, 20 threads: <code>results/main_t20.txt</code></summary>

```
main, Julia threads = 20
main: timestep_floe_properties!               n=100000  median    112.597 ms  min    102.634 ms  (1126.0 ns/floe, 19 samples)
Peak host memory (max RSS): 1.75 GB
```

</details>

<details>
<summary>branch, 1 thread (<code>BENCH_BACKENDS=CPU</code>): <code>results/branch_t1.txt</code></summary>

```
gpu-version branch, Julia threads = 1
CPU vs main, n=1000: all fields isapprox = true, max relative difference = 6.353598990094692e-15
branch CPU Float64: kernels only              n=100000  median     37.590 ms  min     37.219 ms  (375.9 ns/floe, 20 samples)
branch CPU: FixedWidthFloes(floes)            n=100000  median     56.313 ms  min     55.645 ms  (563.1 ns/floe, 20 samples)
branch CPU: adapt to device                   n=100000  median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 20 samples)
branch CPU: adapt to host                     n=100000  median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 20 samples)
branch CPU: update_floes!                     n=100000  median     13.983 ms  min     13.536 ms  (139.8 ns/floe, 20 samples)
branch CPU Float64: full step as in run!      n=100000  median     82.673 ms  min     81.860 ms  (826.7 ns/floe, 20 samples)
Peak host memory (max RSS): 1.42 GB
```

</details>

<details>
<summary>branch, 20 threads: <code>results/branch_t20.txt</code></summary>

```
NVIDIA GeForce RTX 3050 Ti Laptop GPU, CUDA runtime 13.4.0
gpu-version branch, Julia threads = 20
CPU vs main, n=1000: all fields isapprox = true, max relative difference = 6.353598990094692e-15
CUDA vs main, n=1000: all fields isapprox = true, max relative difference = 6.353598990094692e-15
branch CPU Float64: kernels only              n=100000  median      5.186 ms  min      4.499 ms  (51.9 ns/floe, 20 samples)
branch CPU: FixedWidthFloes(floes)            n=100000  median     66.267 ms  min     63.710 ms  (662.7 ns/floe, 20 samples)
branch CPU: adapt to device                   n=100000  median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 20 samples)
branch CPU: adapt to host                     n=100000  median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 20 samples)
branch CPU: update_floes!                     n=100000  median     19.934 ms  min     19.091 ms  (199.3 ns/floe, 20 samples)
branch CPU Float64: full step as in run!      n=100000  median     54.883 ms  min     51.249 ms  (548.8 ns/floe, 20 samples)
branch CUDA Float64: kernels only             n=100000  median      6.683 ms  min      6.594 ms  (66.8 ns/floe, 20 samples)
branch CUDA Float32: kernels only             n=100000  median      1.759 ms  min      1.676 ms  (17.6 ns/floe, 20 samples)
branch CUDA: FixedWidthFloes(floes)           n=100000  median     52.765 ms  min     50.015 ms  (527.6 ns/floe, 20 samples)
branch CUDA: adapt to device                  n=100000  median     12.299 ms  min     12.147 ms  (123.0 ns/floe, 20 samples)
branch CUDA: adapt to host                    n=100000  median     25.499 ms  min     23.975 ms  (255.0 ns/floe, 20 samples)
branch CUDA: update_floes!                    n=100000  median     20.019 ms  min     17.901 ms  (200.2 ns/floe, 20 samples)
branch CUDA Float64: full step as in run!     n=100000  median     88.780 ms  min     80.823 ms  (887.8 ns/floe, 20 samples)
Peak host memory (max RSS): 2.74 GB
Peak CUDA memory in use (all processes, measured after each benchmark): 0.32 GB
```

</details>

## How to reproduce

```sh
REPO=/path/to/Subzero.jl            # checked out at the gpu-version commit
WORK=$(mktemp -d)                   # work directory for environments, input, and results
$REPO/benchmark/timestep_floe_properties/setup.sh "$WORK" ae96ff2   # worktree of main + 2 envs
cd "$WORK" && $REPO/benchmark/timestep_floe_properties/run_all.sh     # about 7 min
cat results/*.txt
git -C $REPO worktree remove "$WORK/subzero-main"                     # clean up afterwards
```

`run_all.sh` does the following, with `BENCH_THREADS` threads (default: all CPUs):

1. `make_floes.jl 10000 100 floes_10k.jld2` generates the input (about 1 min). This is
   skipped if `floes_10k.jld2` exists in the work directory. To compare against these
   results, copy the input they used (`floes_10k.jld2` in this folder, 24 MB, not
   committed) into the work directory first.
2. `bench_main.jl` with the `main` environment, with 1 and `BENCH_THREADS` threads. `main`
   is a baseline that doesn't change when the branch changes, so this is skipped if
   `results/main_t<threads>.txt` exists. Delete those files to rerun it.
3. `bench_branch.jl` with the branch environment, with 1 thread (CPU only,
   `BENCH_BACKENDS=CPU`) and `BENCH_THREADS` threads (CPU and GPU).

To benchmark 1,000,000 floes on an HPC system, run e.g.
`BENCH_SIZES=1000000 BENCH_THREADS=32 run_all.sh` (about 12 GB of memory).

`setup.sh` adds CUDA.jl to the branch environment if `nvidia-smi` is found, and AMDGPU.jl
if `rocminfo` or `rocm-smi` is found. To choose yourself, set `BENCH_GPU`, for example
`BENCH_GPU=AMDGPU setup.sh "$WORK"` (or `BENCH_GPU=""` for no GPU package).

Environment variables for the benchmark scripts:

- `BENCH_SIZES`: comma-separated numbers of floes (default `100000`)
- `BENCH_SAMPLES`: samples per benchmark (default 20)
- `BENCH_SECONDS`: maximum time per benchmark in s (default 60)
- `BENCH_BACKENDS`: comma-separated, from `CPU`, `CUDA`, `AMDGPU` (default `CPU` plus each of
  CUDA.jl and AMDGPU.jl that is in the environment and has a working GPU)
- `BENCH_FLOES`: input file (default `floes_10k.jld2`)

Run `bench_main.jl` before `bench_branch.jl`, because it writes `result_main.jld2` for the
correctness check.
