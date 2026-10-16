# `fracture_floes!`

The fracture check (`60_fracture_check`) and splitting the floes that fracture into 3 pieces,
as in the docs examples (`deform_on = false`). Few floes fracture in this input (about 1 in
2,000), so most of the time is the fracture check. Every sample starts from a `deepcopy` of the
floes (not timed).

No correctness check: the floes are split in a threaded loop that shares the rng, so the
result isn't reproducible.
