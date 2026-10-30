# timestep_collisions!: floe-floe and floe-domain interactions. Every floe is checked against
# every later floe (O(n²)), so the baseline is extrapolated above 100,000 floes.

register_step!(;
    title = "`timestep_collisions!`",
    baseline = function (ctx)
        floes, model, consts, n_init = ghosted_floes(ctx.n)
        settings = CollisionSettings()
        # The interactions of the ghosts are added to their parents, so every sample starts
        # from a deepcopy
        trial = @bench(Subzero.timestep_collisions!(f, $n_init, $(model.domain), $consts,
                $(input_Δt()), $settings, Threads.SpinLock()),
            f = deepcopy($floes))
        record!(ctx, "timestep_collisions!", trial; primary = true)
    end,
    baseline_result = function (n)
        floes, _ = collided_floes(n)
        floes = floes[1:n]  # without the ghosts
        Dict("collision_force" => collect(floes.collision_force),
            "collision_trq" => collect(floes.collision_trq), "num_inters" => collect(floes.num_inters),
            "interactions" => [f.interactions[1:f.num_inters, :] for f in floes])
    end,
    # Which of a floe's ghosts collides first depends on the order the threads get a lock, so
    # the result changes a little (about 1e-7) with the number of threads
    check_rtol = 1e-6,
    extrapolate_above = 100_000,
    extrapolate_from = [25_000, 50_000, 100_000],
)
