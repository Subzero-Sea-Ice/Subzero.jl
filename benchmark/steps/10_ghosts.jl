# add_ghosts!: add copies of the floes that cross a periodic boundary, and remove them again
# after collisions and ridging/rafting, as timestep_sim! does.

function ghosts!(floes, domain)
    n_init_floes = length(floes)
    Subzero.add_ghosts!(floes, domain)
    for i in reverse(n_init_floes+1:length(floes))
        StructArrays.foreachfield(field -> deleteat!(field, i), floes)
    end
    empty!.(floes.ghosts)
    return
end

register_step!(;
    title = "`add_ghosts!` and removing the ghosts",
    # Adding and removing the ghosts leaves the floes as they were, so there's nothing to reset
    baseline = function (ctx)
        floes = load_floes(ctx.n)
        domain = make_model(load_floes(0), ctx.n).domain
        trial = @bench(ghosts!($floes, $domain))
        record!(ctx, "add and remove ghosts", trial; primary = true)
    end,
    baseline_result = function (n)
        floes, _ = ghosted_floes(n)
        Dict("centroid" => collect(floes.centroid), "ghosts" => collect(floes.ghosts))
    end,
)
