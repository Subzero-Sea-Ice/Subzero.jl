# `timestep_sim!` (whole timestep)

One whole timestep of the input simulation, with coupling, fractures (`HiblerYieldCurve`, 3
pieces, no deformation) and ridging/rafting run every timestep (by default, coupling runs
every 10 timesteps and the others are off), so that every step in the overview except welding
is part of it once. Welding is off: welding every timestep merges
most floes within a few timesteps. The time of this step is the denominator of the "share of
timestep" column.

Every sample is the next timestep of the same simulation, as in `run!`, so there is no setup.
The number of floes grows a little between samples because of fractures.

The current code runs with `Simulation(; backend)` on every backend: the ported steps run on
that backend, the others on the CPU. Its time includes copying the floes between host and
device (`95_transfer`).

Correctness check: the number of floes, the total area and mass, and the total momentum after
one timestep of 1,000 floes.

The baseline is extrapolated for more than 100,000 floes, because collisions make it grow
faster than linearly with the number of floes.

## Observations

From the run of 2026-10-15 (100,000 floes, P-cores, input `2f2e3e6e9230`): collisions take
about a third of a timestep, coupling another third and ridging/rafting a fifth. The ported
`timestep_floe_properties!` is less than 1% of the timestep, so the current code isn't
faster yet. The differences between the baseline, current CPU and CUDA times (21–25 s with
12 threads) come mostly from the variation of the collisions (see
[`20_collisions`](20_collisions.md)).
