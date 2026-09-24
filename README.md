![Subzero.jl](https://github.com/Subzero-Sea-Ice/Subzero.jl/blob/main/docs/src/assets/title.gif?raw=true)

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://subzero-sea-ice.github.io/Subzero.jl/dev/)
[![CI](https://github.com/Subzero-Sea-Ice/Subzero.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Subzero-Sea-Ice/Subzero.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![CodeCov](https://codecov.io/gh/Subzero-Sea-Ice/Subzero.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Subzero-Sea-Ice/Subzero.jl)
[![Status](https://www.repostatus.org/badges/latest/active.svg?style=flat-square)](https://www.repostatus.org/#active)


## _Fast and Flexible Sea Ice Dynamics_

Subzero.jl is a native [Julia](https://julialang.org/) discrete-element model (DEM) for exploring fine-scale sea ice dynamics, reimplementing and enhancing MATLAB model [SubZero by Manucharyan and Montemuro](https://doi.org/10.1029/2022MS003247).

- 🚀 Runs over **35 times faster** that original MATLAB model for title simulation!
- 🧩 Modular simulation model makes it easy to **customize simulations**!
  - Enable and disable physical processes such as fracturing, ridging, and welding
  - Choose algorithms for key processes (or add your own!)

## _Documentation_

To learn how to build and run simulations, [check out our documentation and tutorials](https://subzero-sea-ice.github.io/Subzero.jl/dev/)!

## _Installation_

Subzero is a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/) (version 1.10 or later). We recommend using [`JuliaUp`](https://github.com/JuliaLang/juliaup) so it is easy to change versions in the future.

2. Launch Julia and type

    ```julia
    julia> using Pkg

    julia> Pkg.add("Subzero")
    ```

    This installs the latest version that's _compatible with your current environment_.

3. You are ready to use Subzero! Get started with:

    ```julia
    julia> using Subzero
    ```

## _Contributing_

If you’re interested in contributing to the development Subzero, we would love to have you! We welcome all kinds of contributions from bug reports, to documentation, to features, and suggestions. **We can't wait to talk to you.**

Please see [CONTRIBUTING](CONTRIBUTING.md) for more details.

## _Citing_

If you use Subzero.jl as part of your research, teaching, or other activities, we would be grateful if you could cite our work.
We are currently working on a JOSS paper, which will be linked here. If you are ready to publish before that, please reach out to us to discuss citations.

## _Authors_

- Primary Author: [**Skylar Gering (@skygering)**](https://github.com/skygering)

The list of [Subzero contributors](https://github.com/Subzero-Sea-Ice/Subzero.jl/graphs/contributors):

[![Contributers](https://contrib.rocks/image?repo=Subzero-Sea-Ice/Subzero.jl)](https://github.com/Subzero-Sea-Ice/Subzero.jl/graphs/contributors)
