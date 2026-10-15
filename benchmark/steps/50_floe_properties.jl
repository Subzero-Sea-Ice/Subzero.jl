# timestep_floe_properties!: move the floes and update their velocities, stresses and strain
# from the forces of the coupling and the collisions. Ported: KernelAbstractions kernels on
# FixedWidthFloes.

# Fields to compare between the baseline and the current code after one call
const FLOE_PROPERTIES_FIELDS = (:centroid, :height, :mass, :moment, :α, :u, :v, :ξ,
    :p_dxdt, :p_dydt, :p_dudt, :p_dvdt, :p_dξdt, :p_dαdt, :stress_accum, :stress_instant, :strain)

register_step!(;
    title = "`timestep_floe_properties!`",
    float32 = true,
    baseline = function (ctx)
        floes, fs = load_floes(ctx.n), FloeSettings()
        # The function changes nested arrays of the floes in place, so every sample starts
        # from a deepcopy
        trial = @bench(Subzero.timestep_floe_properties!(f, 1, $(input_Δt()), $fs),
            f = deepcopy($floes))
        record!(ctx, "timestep_floe_properties!", trial; primary = true)
    end,
    # Kernels only: the floes are already on the device. Converting them and copying them
    # between host and device is measured in 95_transfer.
    current = function (ctx, backend, FT)
        fs = FloeSettings(FT)
        src = adapt(backend, to_float(FT, Subzero.FixedWidthFloes(load_floes(ctx.n))))
        d = adapt(backend, deepcopy(src))
        trial = @bench(begin
                Subzero.timestep_floe_properties!($d, $(FT(input_Δt())), $fs; backend = $backend)
                KernelAbstractions.synchronize($backend)
            end,
            reset!($d, $src, $backend))
        record!(ctx, "kernels", trial; FT, primary = true)
        free!(d)
        free!(src)
    end,
    baseline_result = function (n)
        floes = load_floes(n)
        Subzero.timestep_floe_properties!(floes, 1, input_Δt(), FloeSettings())
        floe_fields(floes, FLOE_PROPERTIES_FIELDS)
    end,
    # What timestep_sim! does around the call, including the copies between host and device
    current_result = function (n, backend)
        floes, fs = load_floes(n), FloeSettings()
        dev_floes = adapt(backend, Subzero.FixedWidthFloes(floes))
        Subzero.timestep_floe_properties!(dev_floes, input_Δt(), fs; backend)
        KernelAbstractions.synchronize(backend)
        Subzero.update_floes!(floes, adapt(Array, dev_floes))
        floe_fields(floes, FLOE_PROPERTIES_FIELDS)
    end,
)
