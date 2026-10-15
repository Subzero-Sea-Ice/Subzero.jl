# `fracture_floes!`

The fracture check (`60_fracture_check`) and splitting the floes that fracture into 3 pieces,
as in the docs examples (`deform_on = false`). Few floes fracture in this input (about 1 in
2,000), so most of the time is the fracture check. Every sample starts from a `deepcopy` of the
floes (not timed).

`timestep_sim!` now runs the fracture check on the device (see
[`60_fracture_check`](60_fracture_check.md)) and passes the floes that fracture to
`fracture_floes!` with the `frac_idx` keyword, so only the splitting runs on the host. This
step still calls `fracture_floes!` without `frac_idx`, so that it measures the same work as the
baseline: in the current code, that runs the check on the CPU backend, including converting
the floes to `FixedWidthFloes`.

No correctness check: the floes are split in a threaded loop that shares the rng, so the
result isn't reproducible.
