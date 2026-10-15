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

function bench_timestep_sim!(ctx, sim)
    tstep = Ref(0)
    trial = @bench(Subzero.timestep_sim!($sim, t), t = ($tstep[] += 1))
    record!(ctx, "timestep_sim!", trial; primary = true)
end

# Totals after one timestep: the number of floes changes with fractures and ridging
function timestep_sim_result(sim)
    Subzero.timestep_sim!(sim, 1)
    floes = sim.model.floes
    Dict("nfloes" => [length(floes)], "area" => [sum(floes.area)], "mass" => [sum(floes.mass)],
        "momentum" => [sum(floes.mass .* floes.u), sum(floes.mass .* floes.v)])
end

register_step!(;
    title = "`timestep_sim!` (whole timestep)",
    baseline = ctx -> bench_timestep_sim!(ctx,
        input_simulation(ctx.n; timestep_sim_settings(load_floes(ctx.n))...)),
    current = (ctx, backend, FT) -> bench_timestep_sim!(ctx,
        input_simulation(ctx.n; backend, timestep_sim_settings(load_floes(ctx.n))...)),
    baseline_result = n -> timestep_sim_result(input_simulation(n; timestep_sim_settings(load_floes(n))...)),
    current_result = (n, backend) ->
        timestep_sim_result(input_simulation(n; backend, timestep_sim_settings(load_floes(n))...)),
    extrapolate_above = 100_000,
    extrapolate_from = [25_000, 50_000, 100_000],
)
