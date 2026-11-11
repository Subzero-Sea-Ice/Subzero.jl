# timestep_coupling!: forces of the ocean and the atmosphere on the floes, from the floes'
# subfloe points on the grid. One-way coupling is ported: a KernelAbstractions kernel on
# FixedWidthFloes and the model moved to the device with adapt_model.

const COUPLING_FIELDS = (:fxOA, :fyOA, :trqOA, :hflx_factor)

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
    # Kernels only: the floes and the model are already on the device. Converting and
    # copying the floes is measured in 95_transfer. The forces are overwritten every call
    # and nothing else changes, so there's nothing to reset between samples.
    current = function (ctx, backend, FT)
        floes = load_floes(ctx.n)
        model = make_model(floes, ctx.n)
        consts, settings = input_constants(floes), CouplingSettings()
        fs = input_floe_settings(model.grid)
        dev_floes = adapt(backend, Subzero.FixedWidthFloes(floes))
        dev_model = Subzero.adapt_model(model, backend)
        trial = @bench(begin
                Subzero.timestep_coupling!($dev_floes, $dev_model, $model, $(input_Δt()),
                    $consts, $settings, $fs, $backend)
                KernelAbstractions.synchronize($backend)
            end)
        record!(ctx, "kernels", trial; FT, primary = true)
        free!(dev_floes)
    end,
    baseline_result = function (n)
        floes = load_floes(n)
        model = make_model(floes, n)
        Subzero.timestep_coupling!(model, input_Δt(), input_constants(floes), CouplingSettings(),
            input_floe_settings(model.grid))
        floe_fields(model.floes, COUPLING_FIELDS)
    end,
    # What timestep_sim! does around the call, including the copies between host and device
    current_result = function (n, backend)
        floes = load_floes(n)
        model = make_model(floes, n)
        dev_floes = adapt(backend, Subzero.FixedWidthFloes(floes))
        dev_model = Subzero.adapt_model(model, backend)
        Subzero.timestep_coupling!(dev_floes, dev_model, model, input_Δt(), input_constants(floes),
            CouplingSettings(), input_floe_settings(model.grid), backend)
        KernelAbstractions.synchronize(backend)
        Subzero.update_floes!(floes, adapt(Array, dev_floes))
        floe_fields(floes, COUPLING_FIELDS)
    end,
)
