# timestep_welding!: floes that overlap weld together with a probability. Floes are welded
# within bins of the domain: 10 × 10 bins per copy of the floe field, so about 100 floes
# per bin for every number of floes. Every pair of floes in a bin is checked.

weld_settings(n) = WeldSettings(; weld_on = true, Δts = [1], Nxs = [10ntiles(n)], Nys = [10])

function weld!(floes, model, settings, fs)
    Subzero.timestep_welding!(floes, maximum(floes.id), model.grid, model.domain, settings, fs,
        1, input_Δt(), Xoshiro(1))
end

register_step!(;
    title = "`timestep_welding!`",
    baseline = function (ctx)
        floes = load_floes(ctx.n)
        model = make_model(load_floes(0), ctx.n)
        settings, fs = weld_settings(ctx.n), input_floe_settings(model.grid)
        trial = @bench(weld!(f, $model, $settings, $fs), f = deepcopy($floes))
        record!(ctx, "timestep_welding!", trial; primary = true)
    end,
    baseline_result = function (n)
        floes = load_floes(n)
        model = make_model(load_floes(0), n)
        weld!(floes, model, weld_settings(n), input_floe_settings(model.grid))
        floe_fields(floes, (:area, :mass, :centroid))
    end,
)
