# `add_ghosts!` and removing the ghosts

`add_ghosts!` adds a copy of every floe that crosses a periodic boundary, for collisions and
ridging/rafting. After those, `timestep_sim!` removes the ghosts again with `deleteat!` on
every field of the floe `StructArray`. Both are measured together. Adding and removing the
ghosts leaves the floes as they were, so there's no setup.

Correctness check: the centroids and the ghost indices of 1,000 floes after `add_ghosts!`.
