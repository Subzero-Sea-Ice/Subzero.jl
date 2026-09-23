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

- Full suite, the same as CI (about 6 min):
  `julia --project=. -e 'using Pkg; Pkg.test()'`
  The output is long; redirect it to a file and read the `Test Summary` table at the end.
- Single file: test deps (LibGEOS, Test) are in `[extras]` in `Project.toml`, so
  there is no `test/Project.toml`. Use TestEnv.jl, which you install once with
  `julia -e 'using Pkg; Pkg.add("TestEnv")'`:

  ```sh
  julia --project=. -e 'using TestEnv; TestEnv.activate()
      using JLD2, Logging, NCDatasets, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, Test
      import GeometryOps as GO; import GeometryOps.GeoInterface as GI
      include("test/utils.jl"); include("test/test_floe_utils.jl")'
  ```

  The `using` lines are the ones at the top of `test/runtests.jl`.
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

- Uses CUDA.jl and KernelAbstractions. Check `julia --project=. -e 'using CUDA; @show CUDA.functional()'`
  before running GPU code.
