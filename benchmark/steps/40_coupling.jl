# timestep_coupling!: forces of the ocean and the atmosphere on the floes, from the floes'
# subfloe points on the grid.

register_step!(;
    title = "`timestep_coupling!`",
    # Every call overwrites the forces, so there's nothing to reset
    baseline = function (ctx)
        floes = load_floes(ctx.n)
        model = make_model(floes, ctx.n)
        consts, settings = input_constants(floes), CouplingSettings()
        fs = input_floe_settings(model.grid)
        trial = @bench(Subzero.timestep_coupling!($model, $(input_Δt()), $consts, $settings, $fs))
        record!(ctx, "timestep_coupling!", trial; primary = true)
    end,
    baseline_result = function (n)
        floes = load_floes(n)
        model = make_model(floes, n)
        Subzero.timestep_coupling!(model, input_Δt(), input_constants(floes), CouplingSettings(),
            input_floe_settings(model.grid))
        floe_fields(model.floes, (:fxOA, :fyOA, :trqOA, :hflx_factor))
    end,
)
