# timestep_ridging_rafting!: floes that overlap ridge or raft, which moves height and area
# between them and can break them into pieces. Runs after the collisions, on the floes with
# their ghosts.

ridging_settings() = RidgeRaftSettings(; ridge_raft_on = true, Δt = 1)

function ridge_raft!(floes, model, settings, fs)
    max_floe_id = maximum(floes.id)
    pieces_buffer = StructArray{Floe{Float64}}(undef, 0)
    Subzero.timestep_ridging_rafting!(floes, pieces_buffer, model.domain, max_floe_id, settings,
        fs, SimplificationSettings(), input_Δt(), Xoshiro(1))
    return pieces_buffer
end

register_step!(;
    title = "`timestep_ridging_rafting!`",
    baseline = function (ctx)
        floes, model, _ = collided_floes(ctx.n)
        settings, fs = ridging_settings(), input_floe_settings(model.grid)
        trial = @bench(ridge_raft!(f, $model, $settings, $fs), f = deepcopy($floes))
        record!(ctx, "timestep_ridging_rafting!", trial; primary = true)
    end,
    baseline_result = function (n)
        floes, model, _ = collided_floes(n)
        pieces = ridge_raft!(floes, model, ridging_settings(), input_floe_settings(model.grid))
        Dict("height" => collect(floes.height), "area" => collect(floes.area),
            "mass" => collect(floes.mass), "npieces" => [length(pieces)])
    end,
)
