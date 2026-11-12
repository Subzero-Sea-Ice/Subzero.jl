# `add_ghosts!` and removing the ghosts

`add_ghosts!` adds a copy of every floe that crosses a periodic boundary, for collisions and
ridging/rafting. After those, `timestep_sim!` removes the ghosts again with `deleteat!` on
every field of the floe `StructArray`. Both are measured together. Adding and removing the
ghosts leaves the floes as they were, so there's no setup.

Correctness check: the centroids and the ghost indices of 1,000 floes after `add_ghosts!`.

## Ported

`add_ghosts!(floes::FixedWidthFloes, domain, backend; workspace)` and
`remove_ghosts!(floes, n_floes, backend)`. Each of the two passes (east/west, then
north/south) is a kernel that counts the ghosts of each floe, a scan of the counts with
AcceleratedKernels, a copy of the total to the host, and a kernel that makes the ghosts in
the rows the scan gives them. Removing the ghosts lowers `n_floes` and empties the parents'
lists of ghosts with one kernel.

`current` measures the kernels only: the floes are already on the device, and one
`GhostWorkspace` is used for all samples, so nothing is allocated on the device. Every sample
starts from the same floes (reset in the setup), because adding ghosts moves floes whose
centroid is outside of the domain. `add_ghosts!` alone is recorded too.

## Results

`thinkpad-p14s` (Ryzen AI 9 HX PRO 470, Radeon 890M), 100,000 floes, 2026-11-12:

| | add and remove | `add_ghosts!` |
|---|--:|--:|
| baseline, 1 and 8 threads | 21.4 ms | – |
| CPU, 1 thread | 2.03 ms | 1.85 ms |
| CPU, 8 threads | 1.31 ms | 1.20 ms |
| AMDGPU Float64 | 1.96 ms | 1.73 ms |
| AMDGPU Float32 | 1.59 ms | 1.51 ms |

The checks match the baseline exactly on the CPU and the GPU (maximum relative difference 0).

The ported code is 11–16× faster than the baseline on the CPU, but the GPU is no faster
than the CPU. There is little work: only floes at the edges of the domain get ghosts, 805
of 100,000 floes in this input, and 798 of those in the north/south pass, because the floe
field is copied side by side in x. So most of the time is fixed costs and latency:

- Each pass waits for the device twice and copies the total to the host, and launches two
  kernels, a `fill!` and the scan: about 0.2 ms per pass on the GPU.
- On the GPU, most of the north/south pass (about 1.2 ms of it) is `_make_ghosts!`: each of
  the few hundred threads that make ghosts copies whole floes field by field, one memory
  access after another, so the kernel waits for memory, and there are too few threads to
  hide that. Copying with one thread per value instead of per floe would probably help.
- `remove_ghosts!` takes about 0.1–0.3 ms.

The baseline takes as long with 8 threads as with 1, because it isn't threaded.
