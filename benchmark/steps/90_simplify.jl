# simplify_floes!: smooth floes with too many vertices, fuse floes marked for fusion, and
# remove floes marked for removal.

function simplify!(model, fs)
    Subzero.simplify_floes!(model, maximum(model.floes.id), SimplificationSettings(),
        CollisionSettings(), fs, input_Δt(), Xoshiro(1))
end

register_step!(;
    title = "`simplify_floes!`",
    baseline = function (ctx)
        floes = load_floes(ctx.n)
        model = make_model(floes, ctx.n)
        fs = input_floe_settings(model.grid)
        # Every sample starts from a copy of the floes in a new model with the same grid
        new_model(m) = Model(; grid = m.grid, ocean = m.ocean, atmos = m.atmos, domain = m.domain,
            floes = deepcopy(floes))
        trial = @bench(simplify!(m, $fs), m = $new_model($model))
        record!(ctx, "simplify_floes!", trial; primary = true)
    end,
    baseline_result = function (n)
        floes = load_floes(n)
        model = make_model(floes, n)
        simplify!(model, input_floe_settings(model.grid))
        floe_fields(model.floes, (:area, :mass, :centroid))
    end,
)
