# Fracture check (`determine_fractures`)

The first part of `fracture_floes!`: which floes fracture. It updates the criteria
(Hibler's yield curve with the default settings, as in
`FractureSettings(fractures_on = true)`) from the floes' mean height, and then checks for
every floe whether its stress is outside the yield curve. Ported: `determine_fractures!` runs a
kernel on `FixedWidthFloes` and sets `FLAG_FRACTURE`. `timestep_sim!` runs it right after
`timestep_floe_properties!`, before the floes are copied back to the host.

- **baseline**: `determine_fractures(floes, criteria, FloeSettings())`, a single-threaded loop
  that calls `eigvals` and `GO.coveredby` for every floe. It doesn't change the floes, so
  there's no setup.
- **current**: `determine_fractures!(dev_floes, criteria, FloeSettings(FT); backend)` followed
  by `synchronize`, with `dev_floes` already on the device. Like the baseline, this includes
  updating the criteria (the mean height is computed on the device, the polygon on the host)
  and copying the criteria polygon (100 points) to the device. On GPUs also in Float32, with
  `HiblerYieldCurve(Float32; floes)`.

Correctness check: which of 100,000 floes fracture. Few floes fracture in this input (about
1 in 2,000), so the check uses more floes than other steps. The time of the check doesn't
depend on how many floes fracture.

## Observations

From the run of 2026-10-15 (100,000 floes, P-cores, see the overview for the numbers):

- **The kernel is 8× faster than the baseline on 1 CPU thread** (15.5 vs 127 ms), because it
  doesn't allocate, **and 43× faster with 12 threads** (3.3 vs 142 ms): the baseline doesn't
  use threads, and is a little slower in a 12-thread process.
- **On CUDA, Float64 is slower than 12 CPU threads** (6.9 ms): each floe is tested against
  all 100 edges of the yield curve, and this GPU is slow at Float64. In Float32 it takes
  1.3 ms.
