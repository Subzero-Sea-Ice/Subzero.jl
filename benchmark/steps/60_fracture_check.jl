# determine_fractures: which floes fracture, the first part of fracture_floes!. Updates the
# criteria (Hibler's yield curve with the default settings, as in
# FractureSettings(fractures_on = true)) from the floes' mean height and checks every floe.

register_step!(;
    title = "fracture check (`determine_fractures`, part of `fracture_floes!`)",
    # It doesn't change the floes, so there's nothing to reset
    baseline = function (ctx)
        floes, fs = load_floes(ctx.n), FloeSettings()
        criteria = HiblerYieldCurve(floes)
        trial = @bench(Subzero.determine_fractures($floes, $criteria, $fs))
        record!(ctx, "determine_fractures", trial; primary = true)
    end,
    # Few floes fracture in this input (about 1 in 2000), so check more floes than usual
    baseline_result = function (n)
        floes = load_floes(n)
        Dict("frac_idx" => Subzero.determine_fractures(floes, HiblerYieldCurve(floes), FloeSettings()))
    end,
    check_n = 100_000,
)
