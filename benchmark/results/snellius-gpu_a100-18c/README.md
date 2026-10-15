# Benchmark results: `snellius-gpu_a100-18c`

Made by `summarize.jl` from the results files in this folder: don't edit by hand. See
[the benchmark README](../../README.md) for what the columns mean.

| step | ported | baseline 1 thread | baseline 18 threads | current CPU 1 thread | current CPU 18 threads | current CUDA Float64 | current CUDA Float32 | speedup | share of timestep | check | measured |
|---|:-:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|---|
| [`timestep_sim!` (whole timestep)](../../steps/00_timestep_sim.md) | partly | 85,550 | 30,190 | 70,440 | 30,190 | 30,420 | – | 1× | 100% | ✓ | 2026-10-15 `79e1bf1` |
| [`add_ghosts!` and removing the ghosts](../../steps/10_ghosts.md) | no | 25.1 | 26.8 | 24.6 | 25.0 | – | – | – | <1% | ✓ | 2026-10-15 `79e1bf1` |
| [`timestep_collisions!`](../../steps/20_collisions.md) | no | 56,540 | 6,355 | 58,760 | 6,210 | – | – | – | 21% | ✓ | 2026-10-15 `79e1bf1` |
| [`timestep_ridging_rafting!`](../../steps/30_ridging_rafting.md) | no | 5,464 | 7,351 | 5,432 | 7,399 | – | – | – | 24% | ✓ | 2026-10-15 `79e1bf1` |
| [`timestep_coupling!`](../../steps/40_coupling.md) | no | 16,120 | 16,130 | 16,120 | 15,970 | – | – | – | 53% | ✓ | 2026-10-15 `79e1bf1` |
| [`timestep_floe_properties!`](../../steps/50_floe_properties.md) | yes | 532 | 51.4 | 87.6 | 9.25 | 0.38 | 0.33 | 134× | <1% | ✓ | 2026-10-15 `79e1bf1` |
| [fracture check (`determine_fractures`, part of `fracture_floes!`)](../../steps/60_fracture_check.md) | no | 196 | 217 | 212 | 220 | – | – | – | 1% | ✓ | 2026-10-15 `79e1bf1` |
| [`fracture_floes!` (including the fracture check)](../../steps/70_fracture_floes.md) | no | 250 | 266 | 259 | 269 | – | – | – | 1% | – | 2026-10-15 `79e1bf1` |
| [`timestep_welding!`](../../steps/80_welding.md) | no | 1,462 | 1,415 | 1,478 | 1,436 | – | – | – | 5% | ✓ | 2026-10-15 `79e1bf1` |
| [`simplify_floes!`](../../steps/90_simplify.md) | no | 5.00 | 6.68 | 5.11 | 5.01 | – | – | – | <1% | ✓ | 2026-10-15 `79e1bf1` |
| [host ↔ device: `FixedWidthFloes`, `adapt`, `update_floes!`](../../steps/95_transfer.md) | yes | – | – | 97.2 | 99.4 | 121 | – | – | – | – | 2026-10-15 `79e1bf1` |

Median time of one call in ms for 100,000 floes. ≈: extrapolated from fewer floes as t ∝ nᵇ, with the exponent b.

Machine `snellius-gpu_a100-18c`: Intel(R) Xeon(R) Platinum 8360Y CPU @ 2.40GHz, 1 socket, hyperthreading off, run on CPUs 18-35; GPU: NVIDIA A100-SXM4-40GB; input `unknown` (not the reference input `2f2e3e6e9230`, so not directly comparable).
