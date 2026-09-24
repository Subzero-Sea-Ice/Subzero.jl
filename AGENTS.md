# Subzero.jl — agent guide

Subzero.jl is a Julia package for simulating sea-ice floes. The source is in `src/`,
tests are in `test/`, and docs and demos are in `docs/`. The Makie plotting extension
is in `ext/`. Read CONTRIBUTING.md for code conventions: prefer multiple dispatch
over if/else flags.

## Code layout and style

- `src/Subzero.jl` `include`s every source file in dependency order. A new file
  isn't picked up until it's added there, in the right place.
- `docs/src/api.md` lists functions by hand in `@docs` blocks. Add every new
  public function there. The docs build only warns if one is missing.
- There's no formatter configured. Match the style of the surrounding code and
  don't reformat whole files.

## Setup

- Julia ≥ 1.10 (installed via juliaup). Instantiate once:
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Julia precompiles on first use, which can take several minutes (the docs
  environment takes about 3). Use long command timeouts (10 min) and don't
  interrupt it. Prefer one Julia invocation that does several things over many
  short ones.

## Tests

- Full suite, the same as CI (about 4 min):
  `julia --project=. -e 'using Pkg; Pkg.test()'`
  The output is long; redirect it to a file and read the `Test Summary` table at the end.
- Single file: test deps (LibGEOS, Test) are in `[extras]` in `Project.toml`, so
  there is no `test/Project.toml`. Use TestEnv.jl, which you install once with
  `julia -e 'using Pkg; Pkg.add("TestEnv")'`:

  ```sh
  julia --project=. -e 'using TestEnv; TestEnv.activate()
      using JLD2, Logging, NCDatasets, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, Test
      import GeometryOps as GO; import GeometryOps.GeoInterface as GI; import CUDA
      include("test/utils.jl"); include("test/test_floe_utils.jl")'
  ```

  The `using` lines are the ones at the top of `test/runtests.jl`. Some tests read
  files from `test/inputs/` with relative paths, so run those from `test/` with
  `--project=..` and `include("utils.jl")`.
- When adding a test file, also `include` it from `test/runtests.jl`.
- `test/qualitative_behavior.jl` and `test/compare_results.jl` aren't part of the
  test suite. The first runs simulations you check by eye. The second is out of
  date: it needs the MAT package and MATLAB output that aren't in this repo.
  Don't run or fix them unless asked.
- CI tests on Linux and Windows with the latest Julia release only. Windows is
  checked only in CI, and Julia 1.10 (the minimum in `Project.toml`) isn't
  tested at all.

## Demos and docs

- One-time setup of the docs environment (the same as CI):
  `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`
  This adds Subzero to `docs/Project.toml`. Keep that change locally, but don't commit it.
- The demos are Literate scripts in `docs/literate/examples/` (shear_flow,
  simple_strait, forcing_contained_floes, moving_bounds, restart_sim) and
  `docs/literate/tutorial.jl`. Each one runs a simulation and writes `.jld2`
  files and an `.mp4` into a folder under the current directory, so run them
  from a scratch directory rather than the repo root:

  ```sh
  cd "$(mktemp -d)" && julia --project=/path/to/Subzero.jl/docs /path/to/Subzero.jl/docs/literate/examples/shear_flow.jl
  ```

  shear_flow takes about 1 min once the environment is precompiled.
- Build the docs from the repo root (about 2 min):
  `julia --project=docs docs/make.jl`
  This converts the Literate scripts to Markdown without running them.
  The output goes to `docs/build/`. `deploydocs` is skipped outside CI.

## GPU

- `timestep_sim!` is being ported to KernelAbstractions kernels one function at a
  time. So far `timestep_floe_properties!` is done. The backend comes from
  `Simulation(; backend)` (default `CPU()`), and the same kernels run on the CPU.
  Subzero doesn't load CUDA: users `using CUDA` and pass `CUDABackend()`.
- Check `julia --project=. -e 'using CUDA; @show CUDA.functional()'` before running
  GPU code. The tests skip the CUDA backend when it's `false`.

### Porting a function

- Floes go to the device as a `FixedWidthFloes` (`floe.jl`): plain arrays with the
  floe index first. Ragged fields are padded: `poly` has `n_points` vertices per floe
  and `FILL_VALUE` padding, and `interactions` has `num_inters` rows per floe and zero
  padding. Add new fields there. In `update_floes!`, write changed fields back
  through the StructArray columns (`floes.poly[i] = …`). Iterating over `floes` gives
  copies, so changes made to them are lost.
- Put the per-floe logic in a plain function that takes `(floes::FixedWidthFloes, i, …)`,
  and make it a method of the CPU function when there is one (e.g. `calc_strain!`).
  Don't write a `@kernel` for it: run it with `launch_per_floe!(f, backend, floes, args...)`,
  which calls `f(floes, i, args...)` for every floe from one generic kernel. That way the
  function can be unit-tested on the CPU.
- Test against the CPU version: add a unit test per function, and make sure the
  reference test in `test/test_physical_processes/test_update_floe.jl` still
  passes. It compares against a copy of the original CPU code, on `CPU()` and on
  `CUDABackend()`.
- Keep the numerics the same as the CPU version, even when they look wrong, and report
  suspected bugs instead of fixing them in the port. For example, `calc_strain!` uses
  `u` where `v` is expected for `v1` and `v2`. Fixing it would change results and
  break the reference test.

### Inside kernels

- Don't allocate: no comprehensions, array literals, slices (`A[i, :]`),
  broadcasting or `zeros`. Loop over scalars instead.
- Don't call GeometryOps or GeoInterface, because they build geometry objects. Write
  the geometry maths by hand (see `_move_floe!`). GeometryOps is fine on the CPU.
- Don't log. Set a bit in `floes.flags` (the `FLAG_*` constants in
  `update_floe.jl`). The host logs each type of event once per timestep after copying
  back (`_log_flags`).
- Only pass isbits arguments: scalars or small structs such as stress calculators,
  not `FloeSettings` or `Simulation`.
- `InvalidIRError: ... unsupported dynamic function invocation` means something in
  the kernel isn't GPU-compatible: an allocation, a type instability, or a call into
  a library such as GeoInterface. The stack trace shows where.

### Floating point and reproducibility

- `FILL_VALUE` is a Float64. With Float32 floes, compare the padding with `FT(FILL_VALUE)`.
- Float64 literals such as `1.5Δt` turn Float32 calculations into Float64. That's
  kept for now so that results stay within round-off of `main`, but it is slow on GPUs.
  Consumer GPUs are much slower at Float64 than Float32, so check correctness in
  Float64 and measure speed in Float32.
- A Voronoi floe field (`VoronoiTesselationFieldGenerator`) isn't reproducible with
  a fixed `rng`. To compare two runs, generate the floes once and deep-copy them (see
  `test/test_simulation.jl`). In simulations with collisions, differences at
  round-off level grow over time, so compare early timesteps tightly and later ones
  loosely.
