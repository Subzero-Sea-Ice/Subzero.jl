# Benchmark results: `thinkpad-p14s`

Made by `summarize.jl` from the results files in this folder: don't edit by hand. See
[the benchmark README](../../README.md) for what the columns mean.

| step | ported | baseline 1 thread | baseline 8 threads | current CPU 1 thread | current CPU 8 threads | current AMDGPU Float64 | current AMDGPU Float32 | speedup | share of timestep | check | measured |
|---|:-:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|---|
| [`add_ghosts!` and removing the ghosts](../../steps/10_ghosts.md) | yes | 21.4 | 21.4 | 3.74 | 1.78 | 1.28 | 1.06 | 16.7× | – | ✓ | 2026-11-22 `aab9b8f` |

Median time of one call in ms for 100,000 floes. ≈: extrapolated from fewer floes as t ∝ nᵇ, with the exponent b.

Machine `thinkpad-p14s`: AMD Ryzen AI 9 HX PRO 470 w/ Radeon 890M, 1 socket, hyperthreading on, run on CPUs 0-3,12-15; GPU: AMD Radeon Graphics; input `63619f961b19` (not the reference input `2f2e3e6e9230`, so not directly comparable).
