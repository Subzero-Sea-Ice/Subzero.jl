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
