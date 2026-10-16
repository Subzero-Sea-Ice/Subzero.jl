# `timestep_floe_properties!`

Moves the floes and updates their velocities, stresses and strain from the forces of the
coupling and the collisions. Ported: KernelAbstractions kernels on `FixedWidthFloes`.

- **baseline**: `timestep_floe_properties!(floes, 1, Δt, FloeSettings())`, a
  `Threads.@threads` loop over the floe `StructArray`. Every sample starts from a `deepcopy`
  of the floes (not timed).
- **current**: `timestep_floe_properties!(dev_floes, Δt, FloeSettings(FT); backend)` followed
  by `KernelAbstractions.synchronize(backend)`, with `dev_floes` already on the device. Every
  sample starts by copying the input into `dev_floes` on the device (not timed). Converting
  the floes and copying them between host and device is measured in
  [`95_transfer`](95_transfer.md).
- **current, Float32** (GPU only): the same, with every `Float64` array of the
  `FixedWidthFloes` converted to `Float32` and `FloeSettings(Float32)`. The kernels still do
  some `Float64` maths because of `Float64` literals (`1.5Δt`), see `AGENTS.md`.

Correctness check: the centroids, heights, masses, moments, α, velocities, previous-timestep
values, stresses and strain of 1,000 floes after one call, including the copies between host
and device as in `timestep_sim!`.

## Observations

From the run of 2026-10-15 (100,000 floes, P-cores, see the overview for the numbers):

- **The kernels are 5× faster than the baseline** with 12 threads on the CPU backend (9.3 vs
  49 ms), and 7× with 1 thread (39 vs 263 ms), because they don't allocate. With 12 threads,
  the kernels varied between 6 and 9 ms from run to run.
- **CUDA in Float64 is no faster than 12 CPU threads** (6.6 ms): this GPU is slow at Float64.
  In Float32 it's 4× faster (1.6 ms).
- **Copying the floes costs more than the kernels**: the round trip of
  [`95_transfer`](95_transfer.md) takes 78 ms on CUDA. The GPU only pays off once the floes
  stay on the device between steps.
