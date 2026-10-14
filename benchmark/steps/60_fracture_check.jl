# The fracture check: which floes fracture, the first part of fracture_floes!. Updates the
# criteria (Hibler's yield curve with the default settings, as in
# FractureSettings(fractures_on = true)) from the floes' mean height and checks every floe.
# Ported: determine_fractures! runs a kernel on FixedWidthFloes and sets FLAG_FRACTURE.

register_step!(;
    title = "fracture check (`determine_fractures`, part of `fracture_floes!`)",
    float32 = true,
    # It doesn't change the floes (except their flags on the device), so there's nothing to reset
    baseline = function (ctx)
        floes, fs = load_floes(ctx.n), FloeSettings()
        criteria = HiblerYieldCurve(floes)
        trial = @bench(Subzero.determine_fractures($floes, $criteria, $fs))
        record!(ctx, "determine_fractures", trial; primary = true)
    end,
    # The floes are already on the device. This includes updating the criteria (the mean
    # height on the device, the polygon on the host) and copying the criteria polygon to the
    # device, as in timestep_sim!.
    current = function (ctx, backend, FT)
        floes, fs = load_floes(ctx.n), FloeSettings(FT)
        d = adapt(backend, to_float(FT, Subzero.FixedWidthFloes(floes)))
        criteria = HiblerYieldCurve(FT; floes)
        trial = @bench(begin
            Subzero.determine_fractures!($d, $criteria, $fs; backend = $backend)
            KernelAbstractions.synchronize($backend)
        end)
        record!(ctx, "determine_fractures!", trial; FT, primary = true)
        free!(d)
    end,
    # Few floes fracture in this input (about 1 in 2000), so check more floes than usual
    baseline_result = function (n)
        floes = load_floes(n)
        Dict("frac_idx" => Subzero.determine_fractures(floes, HiblerYieldCurve(floes), FloeSettings()))
    end,
    current_result = function (n, backend)
        floes = load_floes(n)
        d = adapt(backend, Subzero.FixedWidthFloes(floes))
        Subzero.determine_fractures!(d, HiblerYieldCurve(floes), FloeSettings(); backend)
        KernelAbstractions.synchronize(backend)
        flags = Array(d.flags)
        free!(d)
        Dict("frac_idx" => findall(f -> f & Subzero.FLAG_FRACTURE != 0, flags))
    end,
    check_n = 100_000,
)
