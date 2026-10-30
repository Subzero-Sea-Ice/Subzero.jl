# Host ↔ device: `FixedWidthFloes`, `adapt`, `update_floes!`

What `timestep_sim!` does around the ported steps: convert the floe `StructArray` to
`FixedWidthFloes`, copy it to the device with `adapt(backend, …)`, copy it back with
`adapt(Array, …)`, and write the fields back into the `StructArray` with `update_floes!`.
This cost goes away once the floes stay on the device between timesteps. The baseline has
nothing like it.

The primary measurement is the round trip. The parts are measured separately too (in the
results file): `FixedWidthFloes(floes)`, adapt to device, adapt to host, `update_floes!`.
Every sample starts by resetting the floes with `update_floes!` (not timed). On the CPU
backend, `adapt` doesn't copy.

## Observations

From the run of 2026-10-15 (100,000 floes, 12 threads on the P-cores):

| backend | round trip | `FixedWidthFloes` | to device | to host | `update_floes!` |
|---|--:|--:|--:|--:|--:|
| CPU  | 46.6 | 30.6 | 0    | 0    | 16.7 |
| CUDA | 78.0 | 31.0 | 11.9 | 23.0 | 17.6 |

Median times in ms; the parts are measured separately. On CUDA, the round trip takes 12× as
long as the `timestep_floe_properties!` kernels (6.6 ms).
