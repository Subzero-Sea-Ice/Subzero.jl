# `timestep_ridging_rafting!`

Floes that overlap enough ridge or raft: height and area move from one floe to the other,
and floes can break into pieces. It runs after the collisions, on the floes with their
ghosts, with `RidgeRaftSettings(ridge_raft_on = true)` and default settings otherwise. Every
sample starts from a `deepcopy` of the floes after `timestep_collisions!` and a new
`Xoshiro(1)` (not timed).

Correctness check: the heights, areas and masses of 1,000 floes, and the number of pieces.

The input was spun up without ridging, so many floes overlap enough to ridge or raft in the
first call. In a simulation that ridges every timestep, later calls have less to do.

## Observations

The baseline doesn't use threads, and is slower in a 12-thread Julia process than in a
1-thread one (4.8 vs 4.1 s on 2026-10-15; 4.6 vs 2.8 s in an earlier run with another input).
The other steps that don't use threads (ghosts, the fracture check, simplification) are
often slower with 12 threads too. This hasn't been investigated.
