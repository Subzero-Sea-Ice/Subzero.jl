# Subzero.jl — agent guide

Subzero.jl is a Julia package for simulating sea-ice floes. The source is in `src/`,
tests are in `test/`, and docs and demos are in `docs/`. The Makie plotting extension
is in `ext/`. Read CONTRIBUTING.md for code conventions: prefer multiple dispatch
over if/else flags.

## Code layout and style

- `src/Subzero.jl` `include`s every source file in dependency order. A new file
  isn't picked up until it's added there, in the right place. The `GeometryOpsGPU`
  submodule includes its own files, in `src/GeometryOpsGPU/GeometryOpsGPU.jl`.
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
- Single file: test deps (JLArrays, LibGEOS, Pkg, Test) are in `[extras]` in
  `Project.toml`, so there is no `test/Project.toml`. Use TestEnv.jl, which you install once with
  `julia -e 'using Pkg; Pkg.add("TestEnv")'`:

  ```sh
  julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/gpu_packages.jl")
      using JLD2, Logging, NCDatasets, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, Test
      import GeometryOps as GO; import GeometryOps.GeoInterface as GI; import JLArrays
      include("test/utils.jl"); include("test/test_floe_utils.jl")'
  ```

  The `include` and `using` lines are the ones at the top of `test/runtests.jl`. Some tests read
  files from `test/inputs/` with relative paths, so run those from `test/` with
  `--project=..` and `include("gpu_packages.jl")` and `include("utils.jl")`.
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
  time. So far `timestep_floe_properties!` is done. Adding and removing ghosts
  (`add_ghosts!` and `remove_ghosts!` for `FixedWidthFloes`, in `collisions.jl`) is
  ported too, but `timestep_sim!` doesn't use it yet. The backend comes from
  `Simulation(; backend)` (default `CPU()`), and the same kernels run on the CPU.
- The goal is to run on any GPU that KernelAbstractions supports (NVIDIA, AMD, Apple,
  Intel), not just CUDA. Subzero doesn't depend on any GPU package: users load one
  (`using CUDA`, `using AMDGPU`, `using Metal`, …) and pass its backend
  (`CUDABackend()`, `ROCBackend()`, `MetalBackend()`, …). Don't add a GPU package to
  `[deps]` or call vendor-specific functions (`CUDA.@cuda`, `CuArray`, …) in `src/`.
  Move data with `adapt(backend, x)` and `adapt(Array, x)`, and wait with
  `KernelAbstractions.synchronize(backend)`. If vendor-specific code is ever needed,
  put it in a package extension in `ext/`.
- Metal and some Intel GPUs don't support Float64 at all, so kernels must be able to
  run on Float32 floes without any Float64 maths (see the note on Float64 literals below).
- The tests run the kernels on every backend in `test_backends()` (`test/utils.jl`):
  `CPU()`, `JLBackend()` from JLArrays (a reference GPU backend that runs on the CPU,
  so CI tests the GPU code path without a GPU), and the GPUs listed in the test argument
  `--gpus=…` (`Pkg.test(test_args = ["--gpus=CUDA"])`) or, if that isn't given, in the
  environment variable `SUBZERO_TEST_GPUS` (`CUDA` and/or `AMDGPU`, comma-separated). The GPU
  packages are in `[extras]` and `[compat]` but not in the test target, so not everyone
  installs them. `test/gpu_packages.jl` `Pkg.add`s the listed ones to the test environment,
  and `test/utils.jl` errors if one can't use a GPU. Include `gpu_packages.jl` before loading
  Subzero: adding a GPU package can change the versions of packages Subzero uses, and a
  loaded package can't be replaced. With current versions they usually line up (AMDGPU 2.9
  and later use LLVM 10, like KernelAbstractions 0.10), but an older one can still clash:
  AMDGPU 2.8 needs LLVM 9, while a resolve without it picks LLVM 10 via KernelAbstractions →
  Atomix → UnsafeAtomics. Check whether a GPU works with
  `SUBZERO_TEST_GPUS=CUDA julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/gpu_packages.jl"); using Subzero; import JLArrays; include("test/utils.jl"); @show test_backends()'`.
  To test another GPU, add its package to `[extras]`, `[compat]` and `GPU_BACKENDS` in
  `test/gpu_packages.jl`. Metal and oneAPI aren't in the list yet: the tests use Float64 floes,
  which Metal and many Intel GPUs can't handle.

### Porting a function

- Floes go to the device as a `FixedWidthFloes` (`floe.jl`): plain arrays with the
  floe index first. Ragged fields are padded: `poly` has `n_points` vertices per floe
  and `interactions` has `num_inters` rows per floe, both padded with zeros. Add new
  fields there. In `update_floes!`, write changed fields back through the StructArray
  columns (`floes.poly[i] = …`). Iterating over `floes` gives copies, so changes made
  to them are lost.
- Put the per-floe logic in a plain function that takes `(floes::FixedWidthFloes, i, …)`,
  and make it a method of the CPU function when there is one (e.g. `calc_strain!`).
  Don't write a `@kernel` for it: run it with `launch_per_floe!(f, backend, floes, args...)`,
  which calls `f(floes, i, args...)` for every floe from one generic kernel. That way the
  function can be unit-tested on the CPU.
- Test against the CPU version: add a unit test per function, and make sure the
  reference test in `test/test_physical_processes/test_update_floe.jl` still
  passes. It compares against a copy of the original CPU code on every backend in
  `test_backends()`.
- Keep the numerics the same as the CPU version, even when they look wrong, and report
  suspected bugs instead of fixing them in the port. For example, `calc_strain!` uses
  `u` where `v` is expected for `v1` and `v2`. Fixing it would change results and
  break the reference test.
- Benchmark it: every step of `timestep_sim!` has a file in `benchmark/steps/`. Add a
  `current` function there, run only that step (`BENCH_STEPS=<step>`), and change only
  that step's files and its row of the overview table (see `benchmark/README.md`).

### Inside kernels

- Don't allocate: no comprehensions, array literals, slices (`A[i, :]`),
  broadcasting or `zeros`. Loop over scalars instead.
- Don't call GeometryOps or GeoInterface, because they build geometry objects. Use the
  `GeometryOpsGPU` submodule (`src/GeometryOpsGPU/`, imported as `GOG`) instead, and add
  missing functions there: e.g. `GOG.transform!(f, _floe_poly(floes, i))` (see
  `_move_floe!`). Its API follows GeometryOps: the same names and argument order where
  GeometryOps has the function, `!` when it writes into its arguments, polygons as
  `GOG.ArrayPolygon`s, and files laid out as in GeometryOps' `src/`. It must not use
  Subzero, so that it can become a package of its own. Test it against GeometryOps in
  `test/test_geometry_ops_gpu.jl`. GeometryOps is fine on the CPU.
- Don't log. Set a bit in `floes.flags` (the `FLAG_*` constants in
  `update_floe.jl`). The host logs each type of event once per timestep after copying
  back (`_log_flags`).
- Only pass isbits arguments: scalars or small structs such as stress calculators,
  not `FloeSettings` or `Simulation`.
- `InvalidIRError: ... unsupported dynamic function invocation` means something in
  the kernel isn't GPU-compatible: an allocation, a type instability, or a call into
  a library such as GeoInterface. The stack trace shows where.
- Inline every function a kernel calls that takes the floes or one of their arrays
  (`@inline` on the method, or at the call). Otherwise they're passed by reference, and on
  CUDA each thread first copies the whole `FixedWidthFloes` into its local memory, which
  makes the kernel several times slower. `test/test_kernel_code.jl` checks this without
  a GPU: add a call that launches each new kernel to it.

### Floating point and reproducibility

- Ported code must match the original (the CPU code on `main`, or the GeometryOps
  function it replaces) bit for bit with Float64 floes, which is how the original was
  run. With Float32 floes it doesn't have to: compare with a tolerance, or not at all.
  This holds for the same function called on the host, without a kernel. In kernels,
  results can differ by round-off: every backend, the CPU backends (POCL) too, uses its own
  `sin`, `cos`, `atan`, ..., which can differ from Julia's by 1–2 ulp, and GPU compilers
  fuse multiplications and additions (CUDA does by default). Which values differ depends
  on the machine, so test kernels that use such functions with a tolerance on every
  backend, also where they match exactly locally. Other kernels (`+`, `*`, `/`, `sqrt`) can
  be tested exactly on the CPU backends.
- With Float32 floes, use Float32 maths throughout. Don't promote to Float64 inside a
  function to copy the original's Float32 results (GeometryOps does this in places, e.g.
  its orientation tests and clipping). Float64 maths is slow on consumer GPUs and fails on
  GPUs without Float64 (Metal), so check correctness in Float64 and measure speed in Float32.
  `test/test_kernel_code.jl` checks that the final LLVM of the kernels has no Float64 with
  Float32 floes.
- Float64 literals such as `1.5Δt` turn Float32 calculations into Float64. Some are
  still in the code; they can be replaced with `FT(1.5) * Δt`, which gives the same
  result with Float64 floes.
- A Voronoi floe field (`VoronoiTesselationFieldGenerator`) isn't reproducible with
  a fixed `rng`. To compare two runs, generate the floes once and deep-copy them (see
  `test/test_simulation.jl`). In simulations with collisions, differences at
  round-off level grow over time, so compare early timesteps tightly and later ones
  loosely.
