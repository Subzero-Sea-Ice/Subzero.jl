# Fracture check (`determine_fractures`)

The first part of `fracture_floes!`: which floes fracture. It updates the criteria
(Hibler's yield curve with the default settings, as in
`FractureSettings(fractures_on = true)`) from the floes' mean height, and then checks for
every floe whether its stress is outside the yield curve. It doesn't change the floes, so
there's no setup.

Correctness check: which of 100,000 floes fracture. Few floes fracture in this input (about
1 in 2,000), so the check uses more floes than other steps. The time of the check doesn't
depend on how many floes fracture.
