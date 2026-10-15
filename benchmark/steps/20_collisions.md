# `timestep_collisions!`

Floe-floe and floe-domain interactions, on the floes with their ghosts. A threaded loop
checks every floe against every later floe (O(n²) bounding-circle tests), and computes the
overlap of the pairs that are close with polygon clipping. A serial pass then copies each
interaction to the other floe, and adds the interactions of the ghosts to their parents.
Every sample starts from a `deepcopy` of the floes with ghosts (not timed).

Correctness check: the collision forces, torques and interactions of 1,000 floes, with a
relative tolerance of 1e-6. The result depends on the number of threads: a floe and its
ghost can both collide with the same floe, and the threaded loop keeps the pair that gets
the lock first. With 1 thread, the current CPU result differed from the 12-thread baseline
by 4e-7 (relative).

The baseline is extrapolated for more than 100,000 floes, because of the O(n²) pair search.

## Observations

From the runs of 2026-10-15 (100,000 floes, P-cores): the time varies a lot between samples
with 12 threads, e.g. a median of 8.7 s and a minimum of 7.0 s, and a median of 12.0 s with a
minimum of 7.5 s in an earlier run. Every floe takes the lock of a shared `Dict` for every
pair that is close. 12 threads are only about 4× faster than 1 thread (37 s).
