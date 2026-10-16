# `simplify_floes!`

Smooths floes with more than 30 vertices, fuses floes marked for fusion, and removes floes
marked for removal (`SimplificationSettings()`). In this input, floes have 7 vertices on
average and none are marked, so this step does little. Every sample starts from a new model
with a `deepcopy` of the floes (not timed).

Correctness check: the areas, masses and centroids of 1,000 floes.
