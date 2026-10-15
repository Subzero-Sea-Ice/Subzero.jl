# `timestep_coupling!`

The forces and torques of the ocean and the atmosphere on the floes, from each floe's
subfloe points on the grid (one-way coupling, `CouplingSettings()`). The grid of the input
has 2 km cells, about 36 per floe, so this step also depends on the size of the grid: for
100,000 floes, 7,080 × 708 cells. Every call overwrites the forces, so there's no setup.
The baseline doesn't use threads: it takes the same time with 1 and 20 threads (see
[CPU cores](../README.md#cpu-cores)).

Correctness check: `fxOA`, `fyOA`, `trqOA` and `hflx_factor` of 1,000 floes.

The input files made before `make_input.jl` saved the subfloe points don't have them. Then
`lib/input.jl` makes them again from the floes' polygons. These points are in the floes'
current orientation instead of their orientation when they were made, which doesn't matter
for the time.
