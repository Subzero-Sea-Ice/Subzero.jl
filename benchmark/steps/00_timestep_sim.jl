# timestep_sim!: one whole timestep, the denominator of the "share" column of the overview.
# Coupling, fractures and ridging/rafting run every timestep (on top of collisions and
# simplification, as in the input simulation), so that every step below is part of it once.
# Welding isn't: welding every timestep merges most floes within a few timesteps.
#
# Every sample is the next timestep of the same simulation, as in run!. The number of floes
# changes a little between samples (fractures).

timestep_sim_settings(floes) = (;
    coupling_settings = CouplingSettings(; Δt = 1),
    fracture_settings = FractureSettings(; fractures_on = true, criteria = HiblerYieldCurve(floes),
        Δt = 1, npieces = 3, deform_on = false),
    ridgeraft_settings = RidgeRaftSettings(; ridge_raft_on = true, Δt = 1),
)

# The simulation of the input, for the baseline and for the current code on `backend`
baseline_sim(n) = input_simulation(n; timestep_sim_settings(load_floes(n))...)
current_sim(n, backend) = input_simulation(n; backend, timestep_sim_settings(load_floes(n))...)

# A function that runs timestep t of `sim`. The current timestep_sim! also takes the grid,
# ocean, atmosphere and domain on the device, made once per simulation as in run!.
baseline_step(sim) = t -> Subzero.timestep_sim!(sim, t)
function current_step(sim)
    dev_model = Subzero.adapt_model(sim.model, sim.backend)
    return t -> Subzero.timestep_sim!(sim, dev_model, t)
end

function bench_timestep_sim!(ctx, step!)
    tstep = Ref(0)
    trial = @bench($step!(t), t = ($tstep[] += 1))
    record!(ctx, "timestep_sim!", trial; primary = true)
end

# Totals after one timestep: the number of floes changes with fractures and ridging
function timestep_sim_result(sim, step!)
    step!(1)
    floes = sim.model.floes
    Dict("nfloes" => [length(floes)], "area" => [sum(floes.area)], "mass" => [sum(floes.mass)],
        "momentum" => [sum(floes.mass .* floes.u), sum(floes.mass .* floes.v)])
end

register_step!(;
    title = "`timestep_sim!` (whole timestep)",
    baseline = ctx -> bench_timestep_sim!(ctx, baseline_step(baseline_sim(ctx.n))),
    current = (ctx, backend, FT) ->
        bench_timestep_sim!(ctx, current_step(current_sim(ctx.n, backend))),
    baseline_result = function (n)
        sim = baseline_sim(n)
        timestep_sim_result(sim, baseline_step(sim))
    end,
    current_result = function (n, backend)
        sim = current_sim(n, backend)
        timestep_sim_result(sim, current_step(sim))
    end,
    extrapolate_above = 100_000,
    extrapolate_from = [25_000, 50_000, 100_000],
)