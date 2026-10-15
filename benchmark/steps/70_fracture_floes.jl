# fracture_floes!: the fracture check (60_fracture_check) and splitting the floes that
# fracture into pieces. Settings as in the docs examples: 3 pieces, no deformation.
#
# No correctness check: the floes are split in a threaded loop that shares the rng, so the
# result isn't reproducible.

fracture_settings(floes) = FractureSettings(; fractures_on = true,
    criteria = HiblerYieldCurve(floes), Δt = 1, npieces = 3, deform_on = false)

register_step!(;
    title = "`fracture_floes!` (including the fracture check)",
    baseline = function (ctx)
        floes = load_floes(ctx.n)
        settings = fracture_settings(floes)
        fs = input_floe_settings(make_model(load_floes(0), ctx.n).grid)
        max_floe_id = maximum(floes.id)
        trial = @bench(Subzero.fracture_floes!(f, $max_floe_id, Xoshiro(1), $settings, $fs,
                $(input_Δt())),
            f = deepcopy($floes))
        record!(ctx, "fracture_floes!", trial; primary = true)
    end,
)
