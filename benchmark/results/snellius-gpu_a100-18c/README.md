# Benchmark results: `snellius-gpu_a100-18c`

Made by `summarize.jl` from the results files in this folder: don't edit by hand. See
[the benchmark README](../../README.md) for what the columns mean.

| step | ported | baseline 1 thread | baseline 18 threads | current CPU 1 thread | current CPU 18 threads | current CUDA Float64 | current CUDA Float32 | speedup | share of timestep | check | measured |
|---|:-:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|---|
| [`timestep_sim!` (whole timestep)](../../steps/00_timestep_sim.md) | partly | 77,540 | 30,380 | 75,090 | 33,400 | 30,880 | – | 0.984× | 100% | ✓ | 2026-10-15 `441c557` |
| [`add_ghosts!` and removing the ghosts](../../steps/10_ghosts.md) | no | 26.2 | 27.6 | 25.5 | 26.1 | – | – | – | <1% | ✓ | 2026-10-15 `441c557` |
| [`timestep_collisions!`](../../steps/20_collisions.md) | no | 53,240 | 7,468 | 64,320 | 7,526 | – | – | – | 25% | ✓ | 2026-10-15 `441c557` |
| [`timestep_ridging_rafting!`](../../steps/30_ridging_rafting.md) | no | 6,562 | 8,108 | 6,608 | 8,046 | – | – | – | 27% | ✓ | 2026-10-15 `441c557` |
| [`timestep_coupling!`](../../steps/40_coupling.md) | no | 16,430 | 15,990 | 16,260 | 15,860 | – | – | – | 53% | ✓ | 2026-10-15 `441c557` |
| [`timestep_floe_properties!`](../../steps/50_floe_properties.md) | yes | 532 | 47.4 | 88.0 | 6.57 | 0.38 | 0.33 | 125× | <1% | ✓ | 2026-10-15 `441c557` |
| [fracture check (`determine_fractures`, part of `fracture_floes!`)](../../steps/60_fracture_check.md) | no | 194 | 213 | 207 | 215 | – | – | – | 1% | ✓ | 2026-10-15 `441c557` |
| [`fracture_floes!` (including the fracture check)](../../steps/70_fracture_floes.md) | no | 292 | 290 | 302 | 282 | – | – | – | 1% | – | 2026-10-15 `441c557` |
| [`timestep_welding!`](../../steps/80_welding.md) | no | 1,474 | 1,418 | 1,485 | 1,456 | – | – | – | 5% | ✓ | 2026-10-15 `441c557` |
| [`simplify_floes!`](../../steps/90_simplify.md) | no | 4.94 | 6.10 | 4.93 | 4.91 | – | – | – | <1% | ✓ | 2026-10-15 `441c557` |
| [host ↔ device: `FixedWidthFloes`, `adapt`, `update_floes!`](../../steps/95_transfer.md) | yes | – | – | 98.6 | 98.9 | 123 | – | – | – | – | 2026-10-15 `441c557` |

Median time of one call in ms for 100,000 floes. ≈: extrapolated from fewer floes as t ∝ nᵇ, with the exponent b.

Machine `snellius-gpu_a100-18c`: Intel(R) Xeon(R) Platinum 8360Y CPU @ 2.40GHz, 1 socket, hyperthreading off, run on CPUs 0-17; GPU: NVIDIA A100-SXM4-40GB; the reference input `2f2e3e6e9230`.
