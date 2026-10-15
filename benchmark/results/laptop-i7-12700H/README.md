# Benchmark results: `laptop-i7-12700H`

Made by `summarize.jl` from the results files in this folder: don't edit by hand. See
[the benchmark README](../../README.md) for what the columns mean.

| step | ported | baseline 1 thread | baseline 12 threads | current CPU 1 thread | current CPU 12 threads | current CUDA Float64 | current CUDA Float32 | speedup | share of timestep | check | measured |
|---|:-:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|---|
| [`timestep_sim!` (whole timestep)](../../steps/00_timestep_sim.md) | partly | 42,990 | 24,750 | 49,860 | 21,010 | 22,520 | – | 1.18× | 100% | ✓ | 2026-10-15 `9d4fe06` |
| [`add_ghosts!` and removing the ghosts](../../steps/10_ghosts.md) | no | 12.1 | 15.1 | 11.6 | 16.5 | – | – | – | <1% | ✓ | 2026-10-15 `9d4fe06` |
| [`timestep_collisions!`](../../steps/20_collisions.md) | no | 37,130 | 8,716 | 43,080 | 8,672 | – | – | – | 35% | ✓ | 2026-10-15 `9d4fe06` |
| [`timestep_ridging_rafting!`](../../steps/30_ridging_rafting.md) | no | 4,054 | 4,780 | 3,925 | 5,061 | – | – | – | 19% | ✓ | 2026-10-15 `9d4fe06` |
| [`timestep_coupling!`](../../steps/40_coupling.md) | no | 7,736 | 8,255 | 7,996 | 8,374 | – | – | – | 33% | ✓ | 2026-10-15 `9d4fe06` |
| [`timestep_floe_properties!`](../../steps/50_floe_properties.md) | yes | 263 | 48.5 | 39.3 | 9.34 | 6.57 | 1.64 | 7.39× | <1% | ✓ | 2026-10-15 `9d4fe06` |
| [fracture check (`determine_fractures`, part of `fracture_floes!`)](../../steps/60_fracture_check.md) | no | 126 | 138 | 129 | 140 | – | – | – | 1% | ✓ | 2026-10-15 `9d4fe06` |
| [`fracture_floes!` (including the fracture check)](../../steps/70_fracture_floes.md) | no | 191 | 175 | 187 | 176 | – | – | – | 1% | – | 2026-10-15 `9d4fe06` |
| [`timestep_welding!`](../../steps/80_welding.md) | no | 791 | 702 | 821 | 680 | – | – | – | 3% | ✓ | 2026-10-15 `9d4fe06` |
| [`simplify_floes!`](../../steps/90_simplify.md) | no | 2.65 | 2.60 | 2.54 | 2.76 | – | – | – | <1% | ✓ | 2026-10-15 `9d4fe06` |
| [host ↔ device: `FixedWidthFloes`, `adapt`, `update_floes!`](../../steps/95_transfer.md) | yes | – | – | 30.8 | 46.6 | 78.0 | – | – | – | – | 2026-10-15 `9d4fe06` |

Median time of one call in ms for 100,000 floes. ≈: extrapolated from fewer floes as t ∝ nᵇ, with the exponent b.

Machine `laptop-i7-12700H`: 12th Gen Intel(R) Core(TM) i7-12700H, 1 socket, hyperthreading on, run on CPUs 0-11; GPU: NVIDIA GeForce RTX 3050 Ti Laptop GPU; the reference input `2f2e3e6e9230`.
