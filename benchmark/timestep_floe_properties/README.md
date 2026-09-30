# Benchmark: `timestep_floe_properties!` on CUDA vs `main`

Compares `timestep_floe_properties!` on the `gpu-version` branch (KernelAbstractions
kernels on `FixedWidthFloes`) with the same function on `main` (a `Threads.@threads`
loop over the floe `StructArray`).

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

Sizes benchmarked: 1,000 floes (the first 1,000), 10,000, and 100,000. For 100,000, the
10,029-floe field is copied 10 times, each copy shifted by `Lx` in x. The function works
on each floe independently, so this only scales the amount of work.

The floe settings are `FloeSettings()` (the default `DecayAreaScaledCalculator`,
`max_floe_height = 10`, `maximum_ξ = 1e-5`, as in the example).

## What is measured

All numbers come from BenchmarkTools (`evals = 1`, up to 15 s per benchmark). Every sample
starts from a fresh `deepcopy` of the input, made in `setup` and not timed.

- **main**: `timestep_floe_properties!(floes, 1, Δt, floe_settings)`. Its per-event
  `@info` logging is silenced with a `NullLogger`.
- **branch, kernels only**: `timestep_floe_properties!(dev_floes, Δt, floe_settings; backend)`
  followed by `KernelAbstractions.synchronize(backend)`, with `dev_floes` already on the
  device. This is 7 kernel launches.
- **branch, Float32 kernels only** (CUDA only): the same, with every `Float64` array of the
  `FixedWidthFloes` converted to `Float32` and `FloeSettings(Float32)`. The kernels still do
  some `Float64` maths because of `Float64` literals (`1.5Δt`), see `AGENTS.md`.
- **branch, full step**: what `run!` does around the call: `FixedWidthFloes(floes)`,
  `adapt(backend, …)`, the kernels, `synchronize`, `adapt(Array, …)`, `update_floes!`,
  and `_log_flags`. It is also measured in parts.

Correctness check: `bench_main.jl` saves the result of one call on 1,000 floes, and
`bench_branch.jl` compares the branch's result with it (centroid, height, mass, moment, α,
velocities, previous-timestep values, stresses, strain). Both CPU and CUDA match `main`,
with a maximum relative difference of 6e-15.

## Machine and versions

- Laptop, on AC power, `balanced` platform profile, `powersave` CPU governor
- CPU: Intel Core i7-12700H (6 P-cores + 8 E-cores, 20 threads)
- GPU: NVIDIA GeForce RTX 3050 Ti Laptop GPU, 4 GB, driver 580.178.04. Consumer GPUs have
  low Float64 throughput (1/64 of Float32 for this one).
- Julia 1.12.7, CUDA.jl 6.4.0 (CUDA runtime 13.4), KernelAbstractions 0.9.43,
  BenchmarkTools 1.8.0
- Subzero: `main` at 4b949fe, `gpu-version` at cba2dcf
- Run on 2026-09-30

## Results

Median time per call in ms (Float64 unless noted). The full output, with minimum times and
sample counts, is under [Example output](#example-output).

| Floes | main, 1 thread | main, 20 threads | branch CPU, 1 thread, kernels | branch CPU, 20 threads, kernels | branch CUDA, kernels | branch CUDA Float32, kernels |
|---:|---:|---:|---:|---:|---:|---:|
| 1,000   |   2.11 |  0.43 |  0.42 | 0.54 | 0.51 | 0.15 |
| 10,000  |  20.7  |  3.09 |  3.85 | 1.01 | 0.97 | 0.62 |
| 100,000 | 364    | 98.4  | 38.7  | 6.96 | 7.81 | 5.01 |

The full step as in `run!`, with the time spent in each part (branch, 20 threads):

| Floes | backend | full step | `FixedWidthFloes` | to device | kernels | to host | `update_floes!` |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1,000   | CPU  |  0.84 |  0.17 | 0     | 0.54 | 0     |  0.13 |
| 1,000   | CUDA |  1.71 |  0.17 | 0.13  | 0.51 | 0.33  |  0.14 |
| 10,000  | CPU  |  4.42 |  2.74 | 0     | 1.01 | 0     |  1.36 |
| 10,000  | CUDA |  7.94 |  2.74 | 1.13  | 0.97 | 1.61  |  1.42 |
| 100,000 | CPU  |  64.5 |  36.7 | 0     | 6.96 | 0     |  15.9 |
| 100,000 | CUDA |   137 |  44.8 | 12.2  | 7.81 | 20.7  |  15.1 |

The 100,000-floe full-step and `update_floes!` numbers are based on only 5–6 samples,
because the untimed `deepcopy` in `setup` uses up most of the time budget. They are noisy:
the CUDA full step (137 ms) is well above the sum of its parts (about 100 ms).

## Observations

- **The kernels are fast.** With 10,000 floes, the CUDA kernels take 0.97 ms, against
  3.1 ms on `main` with 20 threads (3×) and 20.7 ms with 1 thread (21×). With 100,000
  floes it's 7.8 ms against 98 ms (13×) and 364 ms (47×). With 1,000 floes there isn't
  enough work: the time is dominated by 7 kernel launches plus a synchronize (about
  0.5 ms), and `main` is slightly faster.
- **The same kernels on the CPU backend are also much faster than `main`.** On this
  laptop, the kernels with 20 threads are as fast as the Float64 CUDA kernels. With 1
  thread they are 5–9× faster than `main` with 1 thread, because they don't allocate
  (`main` allocates in `calc_stress!`, `maximum(abs.(cforce))`, `_move_floe!`, …).
- **Float64 is slow on this GPU.** Converting the floes to Float32 makes the CUDA kernels
  1.5–3.5× faster, even with the remaining Float64 literals. A datacenter GPU (A100/H100)
  with full-rate Float64 would perform very differently.
- **In `run!`, the conversion costs more than it saves.** Converting the
  `StructArray{Floe}` to `FixedWidthFloes` and back (`update_floes!`) takes 4–8× as long
  as the kernels with 10,000 and 100,000 floes. On CUDA, the host↔device copies (about
  270–330 ns per floe) come on top of that. So the CUDA full step (7.9 ms at 10,000 floes) is currently slower than `main`
  with 20 threads (3.1 ms), and the CPU backend full step (4.4 ms) is too. The GPU only
  pays off once more of `timestep_sim!` runs on the device and the floes stay there
  between steps.
- `main` with 20 threads at 100,000 floes (98 ms, 3.7× faster than 1 thread) scales worse
  than at 10,000 floes (6.7×), probably because of GC pressure from its allocations.

## Example output

The output of `run_all.sh` (the files in `results/`) for the results above.

<details>
<summary><code>main</code>, 1 thread: <code>results/main_t1.txt</code></summary>

```
main, Julia threads = 1
main: timestep_floe_properties!               n=1000    median      2.109 ms  min      1.732 ms  (2108.9 ns/floe, 2104 samples)
main: timestep_floe_properties!               n=10000   median     20.688 ms  min     19.888 ms  (2068.8 ns/floe, 143 samples)
main: timestep_floe_properties!               n=100000  median    364.100 ms  min    254.504 ms  (3641.0 ns/floe, 5 samples)
```

</details>

<details>
<summary><code>main</code>, 20 threads: <code>results/main_t20.txt</code></summary>

```
main, Julia threads = 20
main: timestep_floe_properties!               n=1000    median      0.432 ms  min      0.310 ms  (432.0 ns/floe, 2822 samples)
main: timestep_floe_properties!               n=10000   median      3.089 ms  min      2.157 ms  (308.9 ns/floe, 184 samples)
main: timestep_floe_properties!               n=100000  median     98.423 ms  min     87.357 ms  (984.2 ns/floe, 5 samples)
```

</details>

<details>
<summary>branch, 1 thread (<code>BENCH_BACKENDS=CPU</code>): <code>results/branch_t1.txt</code></summary>

```
NVIDIA GeForce RTX 3050 Ti Laptop GPU, CUDA runtime 13.4.0
gpu-version branch, Julia threads = 1
CPU vs main, n=1000: all fields isapprox = true, max relative difference = 6.353598990094692e-15
branch CPU Float64: kernels only              n=1000    median      0.419 ms  min      0.335 ms  (419.2 ns/floe, 10000 samples)
branch CPU: FixedWidthFloes(floes)            n=1000    median      0.156 ms  min      0.138 ms  (156.2 ns/floe, 10000 samples)
branch CPU: adapt to device                   n=1000    median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: adapt to host                     n=1000    median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: update_floes!                     n=1000    median      0.123 ms  min      0.107 ms  (123.1 ns/floe, 3370 samples)
branch CPU Float64: full step as in run!      n=1000    median      0.784 ms  min      0.654 ms  (783.9 ns/floe, 2789 samples)
branch CPU Float64: kernels only              n=10000   median      3.845 ms  min      3.400 ms  (384.5 ns/floe, 2191 samples)
branch CPU: FixedWidthFloes(floes)            n=10000   median      2.630 ms  min      2.353 ms  (263.0 ns/floe, 3216 samples)
branch CPU: adapt to device                   n=10000   median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: adapt to host                     n=10000   median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: update_floes!                     n=10000   median      1.385 ms  min      1.224 ms  (138.5 ns/floe, 198 samples)
branch CPU Float64: full step as in run!      n=10000   median      7.182 ms  min      6.466 ms  (718.2 ns/floe, 170 samples)
branch CPU Float64: kernels only              n=100000  median     38.699 ms  min     36.907 ms  (387.0 ns/floe, 170 samples)
branch CPU: FixedWidthFloes(floes)            n=100000  median     36.826 ms  min     32.138 ms  (368.3 ns/floe, 217 samples)
branch CPU: adapt to device                   n=100000  median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: adapt to host                     n=100000  median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: update_floes!                     n=100000  median     18.838 ms  min     14.663 ms  (188.4 ns/floe, 6 samples)
branch CPU Float64: full step as in run!      n=100000  median     92.719 ms  min     81.915 ms  (927.2 ns/floe, 5 samples)
```

</details>

<details>
<summary>branch, 20 threads: <code>results/branch_t20.txt</code></summary>

```
NVIDIA GeForce RTX 3050 Ti Laptop GPU, CUDA runtime 13.4.0
gpu-version branch, Julia threads = 20
CPU vs main, n=1000: all fields isapprox = true, max relative difference = 6.353598990094692e-15
CUDA vs main, n=1000: all fields isapprox = true, max relative difference = 6.353598990094692e-15
branch CPU Float64: kernels only              n=1000    median      0.541 ms  min      0.299 ms  (540.8 ns/floe, 10000 samples)
branch CPU: FixedWidthFloes(floes)            n=1000    median      0.167 ms  min      0.141 ms  (166.6 ns/floe, 10000 samples)
branch CPU: adapt to device                   n=1000    median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: adapt to host                     n=1000    median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: update_floes!                     n=1000    median      0.131 ms  min      0.108 ms  (131.3 ns/floe, 3525 samples)
branch CPU Float64: full step as in run!      n=1000    median      0.835 ms  min      0.658 ms  (834.9 ns/floe, 2819 samples)
branch CUDA Float64: kernels only             n=1000    median      0.506 ms  min      0.442 ms  (505.9 ns/floe, 10000 samples)
branch CUDA Float32: kernels only             n=1000    median      0.146 ms  min      0.133 ms  (146.3 ns/floe, 10000 samples)
branch CUDA: FixedWidthFloes(floes)           n=1000    median      0.165 ms  min      0.138 ms  (164.9 ns/floe, 10000 samples)
branch CUDA: adapt to device                  n=1000    median      0.129 ms  min      0.116 ms  (128.5 ns/floe, 10000 samples)
branch CUDA: adapt to host                    n=1000    median      0.327 ms  min      0.297 ms  (326.5 ns/floe, 10000 samples)
branch CUDA: update_floes!                    n=1000    median      0.144 ms  min      0.111 ms  (144.3 ns/floe, 3303 samples)
branch CUDA Float64: full step as in run!     n=1000    median      1.712 ms  min      1.535 ms  (1711.9 ns/floe, 2445 samples)
branch CPU Float64: kernels only              n=10000   median      1.008 ms  min      0.602 ms  (100.8 ns/floe, 4227 samples)
branch CPU: FixedWidthFloes(floes)            n=10000   median      2.739 ms  min      2.427 ms  (273.9 ns/floe, 4250 samples)
branch CPU: adapt to device                   n=10000   median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: adapt to host                     n=10000   median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: update_floes!                     n=10000   median      1.361 ms  min      1.257 ms  (136.1 ns/floe, 213 samples)
branch CPU Float64: full step as in run!      n=10000   median      4.416 ms  min      3.826 ms  (441.6 ns/floe, 186 samples)
branch CUDA Float64: kernels only             n=10000   median      0.973 ms  min      0.847 ms  (97.3 ns/floe, 3888 samples)
branch CUDA Float32: kernels only             n=10000   median      0.620 ms  min      0.528 ms  (62.0 ns/floe, 7462 samples)
branch CUDA: FixedWidthFloes(floes)           n=10000   median      2.736 ms  min      2.392 ms  (273.6 ns/floe, 4141 samples)
branch CUDA: adapt to device                  n=10000   median      1.130 ms  min      1.096 ms  (113.0 ns/floe, 10000 samples)
branch CUDA: adapt to host                    n=10000   median      1.606 ms  min      1.521 ms  (160.6 ns/floe, 6031 samples)
branch CUDA: update_floes!                    n=10000   median      1.416 ms  min      1.285 ms  (141.6 ns/floe, 196 samples)
branch CUDA Float64: full step as in run!     n=10000   median      7.935 ms  min      7.275 ms  (793.5 ns/floe, 172 samples)
branch CPU Float64: kernels only              n=100000  median      6.964 ms  min      4.708 ms  (69.6 ns/floe, 341 samples)
branch CPU: FixedWidthFloes(floes)            n=100000  median     36.698 ms  min     33.172 ms  (367.0 ns/floe, 321 samples)
branch CPU: adapt to device                   n=100000  median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: adapt to host                     n=100000  median      0.000 ms  min      0.000 ms  (0.0 ns/floe, 10000 samples)
branch CPU: update_floes!                     n=100000  median     15.873 ms  min     14.798 ms  (158.7 ns/floe, 5 samples)
branch CPU Float64: full step as in run!      n=100000  median     64.512 ms  min     46.515 ms  (645.1 ns/floe, 5 samples)
branch CUDA Float64: kernels only             n=100000  median      7.807 ms  min      7.680 ms  (78.1 ns/floe, 327 samples)
branch CUDA Float32: kernels only             n=100000  median      5.006 ms  min      4.879 ms  (50.1 ns/floe, 773 samples)
branch CUDA: FixedWidthFloes(floes)           n=100000  median     44.838 ms  min     32.774 ms  (448.4 ns/floe, 295 samples)
branch CUDA: adapt to device                  n=100000  median     12.171 ms  min     11.964 ms  (121.7 ns/floe, 1220 samples)
branch CUDA: adapt to host                    n=100000  median     20.686 ms  min     13.609 ms  (206.9 ns/floe, 558 samples)
branch CUDA: update_floes!                    n=100000  median     15.126 ms  min     14.446 ms  (151.3 ns/floe, 5 samples)
branch CUDA Float64: full step as in run!     n=100000  median    137.070 ms  min    113.584 ms  (1370.7 ns/floe, 5 samples)
```

</details>

## How to reproduce

```sh
REPO=/path/to/Subzero.jl            # checked out at the gpu-version commit
WORK=$(mktemp -d)                   # work directory for environments, input, and results
$REPO/benchmark/timestep_floe_properties/setup.sh "$WORK" 4b949fe   # worktree of main + 2 envs
cd "$WORK" && $REPO/benchmark/timestep_floe_properties/run_all.sh     # about 25 min
cat results/*.txt
git -C $REPO worktree remove "$WORK/subzero-main"                     # clean up afterwards
```

`run_all.sh` does the following:

1. `julia -t 8 --project=env-branch make_floes.jl 10000 100 floes_10k.jld2` generates the
   input (about 1 min). This is skipped if `floes_10k.jld2` exists in the work directory.
   To compare against these results, copy the input they used (`floes_10k.jld2` in this
   folder, 24 MB, not committed) into the work directory first.
2. `bench_main.jl` with the `main` environment, with 1 and 20 threads.
3. `bench_branch.jl` with the branch environment, with 1 thread (CPU only,
   `BENCH_BACKENDS=CPU`) and 20 threads (CPU and CUDA).

Environment variables for the benchmark scripts: `BENCH_SIZES` (comma-separated numbers of
floes, default `1000,10000,100000`), `BENCH_SECONDS` (time budget per benchmark, default
10; `run_all.sh` uses 15), `BENCH_BACKENDS` (default `CPU,CUDA`), and `BENCH_FLOES` (input
file, default `floes_10k.jld2`). Run `bench_main.jl` before `bench_branch.jl`, because it
writes `result_main.jld2` for the correctness check.
