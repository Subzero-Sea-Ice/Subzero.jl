# `timestep_welding!`

Floes that overlap weld together with a probability. The floes are sorted into bins by their
centroid, and every pair of floes in a bin is checked (serially). The bins are 10 × 10 per copy
of the input floe field, so there are about 100 floes per bin for every number of floes:
`WeldSettings(weld_on = true, Nxs = [10 * copies], Nys = [10])`. With fewer, larger bins,
the time grows with the square of the number of floes per bin. Every sample starts from a
`deepcopy` of the floes and a new `Xoshiro(1)` (not timed).

Correctness check: the areas, masses and centroids of 1,000 floes.
