# `add_ghosts!` and removing the ghosts

`add_ghosts!` adds a copy of every floe that crosses a periodic boundary, for collisions and
ridging/rafting. After those, `timestep_sim!` removes the ghosts again with `deleteat!` on
every field of the floe `StructArray`. Both are measured together. Adding and removing the
ghosts leaves the floes as they were, so there's no setup.

Correctness check: the centroids and the ghost indices of 1,000 floes after `add_ghosts!`.

## Ported

`add_ghosts!(floes::FixedWidthFloes, domain, backend; workspace)` and
`remove_ghosts!(floes, n_floes, backend)`. Each of the two passes (east/west, then
north/south) is:

1. a kernel that counts the ghosts of each floe,
2. a scan of the counts with AcceleratedKernels, which gives each floe the rows of its
   ghosts, and a copy of the total to the host,
3. a kernel that plans the ghosts of each floe: source row, ghost id, and shift,
4. a kernel with one thread per value of each ghost, which copies the value from the source
   row (and shifts the centroid and the polygon),
5. a kernel that adds the ghosts to the lists of their parents and swaps parents whose
   centroid is outside of the domain with their last ghost.

Removing the ghosts lowers `n_floes` and empties the parents' lists of ghosts with one
kernel.

`current` measures the kernels only: the floes are already on the device, and one
`GhostWorkspace` is used for all samples, so nothing is allocated on the device. Every sample
starts from the same floes (reset in the setup), because adding ghosts moves floes whose
centroid is outside of the domain. `add_ghosts!` alone is recorded too.

## Results

`thinkpad-p14s` (Ryzen AI 9 HX PRO 470, Radeon 890M), 100,000 floes, 2026-11-13, `4f83ec2`:

| | add and remove | `add_ghosts!` |
|---|--:|--:|
| baseline, 1 and 8 threads | 21.4 ms | – |
| CPU, 1 thread | 4.08 ms | 3.96 ms |
| CPU, 8 threads | 1.85 ms | 1.75 ms |
| AMDGPU Float64 | 1.27 ms | 1.03 ms |
| AMDGPU Float32 | 1.08 ms | 1.17 ms |

The checks match the baseline exactly on the CPU and the GPU (maximum relative difference 0).

There is little work: only floes at the edges of the domain get ghosts, 805 of 100,000
floes in this input, and 798 of those in the north/south pass, because the floe field is
copied side by side in x. So a large part of the time is fixed costs: each pass launches
five kernels, a `fill!`, and the scan, and waits for the device twice. On this GPU, a kernel
launch costs tens of µs.

The ghosts are copied with one thread per value, rather than one thread per floe that
copies a whole floe, value after value. Measured on 2026-11-12, on the GPU, that took the north/south pass from
about 1.2 ms to 0.7 ms, and adding and removing the ghosts from 1.96 ms to 1.31 ms
(Float64). On the CPU it's slower (1.31 ms to 1.93 ms with 8 threads, 2.03 ms to 4.33 ms
with 1 thread), because every thread first has to find the field of its value. The GPU
matters more, and both are much faster than the baseline.

The baseline takes as long with 8 threads as with 1, because it isn't threaded.
