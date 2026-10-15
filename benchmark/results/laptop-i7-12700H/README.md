# Benchmark results: `laptop-i7-12700H`

Made by `summarize.jl` from the results files in this folder: don't edit by hand. See
[the benchmark README](../../README.md) for what the columns mean.

| step | ported | baseline 1 thread | baseline 12 threads | current CPU 1 thread | current CPU 12 threads | current CUDA Float64 | current CUDA Float32 | speedup | share of timestep | check | measured |
|---|:-:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|---|
| [`timestep_sim!` (whole timestep)](../../steps/00_timestep_sim.md) | partly | 45,260 | 17,970 | 54,580 | 18,510 | 19,940 | – | 0.971× | 100% | ✓ | 2026-10-15 `cb4f860` |
| [`add_ghosts!` and removing the ghosts](../../steps/10_ghosts.md) | no | 12.6 | 16.0 | 11.9 | 16.3 | – | – | – | <1% | ✓ | 2026-10-15 `1160b6a` |
| [`timestep_collisions!`](../../steps/20_collisions.md) | no | 35,900 | 9,551 | 34,480 | 12,030 | – | – | – | 53% | ✓ | 2026-10-15 `1160b6a` |
| [`timestep_ridging_rafting!`](../../steps/30_ridging_rafting.md) | no | 2,845 | 4,554 | 2,726 | 4,571 | – | – | – | 25% | ✓ | 2026-10-15 `1160b6a` |
| [`timestep_coupling!`](../../steps/40_coupling.md) | no | 6,156 | 5,399 | 5,488 | 5,211 | – | – | – | 30% | ✓ | 2026-10-15 `1160b6a` |
| [`timestep_floe_properties!`](../../steps/50_floe_properties.md) | yes | 383 | 49.0 | 42.3 | 7.19 | 6.69 | 1.64 | 7.33× | <1% | ✓ | 2026-10-15 `cb4f860` |
| [fracture check (`determine_fractures`, part of `fracture_floes!`)](../../steps/60_fracture_check.md) | yes | 127 | 142 | 15.5 | 3.26 | 6.93 | 1.30 | 43.4× | 1% | ✓ | 2026-10-15 `cb4f860` |
| [`fracture_floes!` (including the fracture check)](../../steps/70_fracture_floes.md) | no | 336 | 178 | 237 | 106 | – | – | – | 1% | – | 2026-10-15 `cb4f860` |
| [`timestep_welding!`](../../steps/80_welding.md) | no | 662 | 705 | 660 | 701 | – | – | – | 4% | ✓ | 2026-10-15 `1160b6a` |
| [`simplify_floes!`](../../steps/90_simplify.md) | no | 2.48 | 2.60 | 2.48 | 2.46 | – | – | – | <1% | ✓ | 2026-10-15 `1160b6a` |
| [host ↔ device: `FixedWidthFloes`, `adapt`, `update_floes!`](../../steps/95_transfer.md) | yes | – | – | 31.3 | 34.5 | 63.1 | – | – | – | – | 2026-10-15 `1160b6a` |

Median time of one call in ms for 100,000 floes. ≈: extrapolated from fewer floes as t ∝ nᵇ, with the exponent b.

Machine `laptop-i7-12700H`: 12th Gen Intel(R) Core(TM) i7-12700H, 1 socket, hyperthreading on, run on CPUs 0-11; GPU: NVIDIA GeForce RTX 3050 Ti Laptop GPU.
