# add_ghosts!: add copies of the floes that cross a periodic boundary, and remove them again
# after collisions and ridging/rafting, as timestep_sim! does. Ported: KernelAbstractions
# kernels on FixedWidthFloes, with a scan from AcceleratedKernels.

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
    #= Kernels only: the floes are already on the device. Every sample starts from the same
    floes, because adding ghosts moves floes whose centroid is outside of the domain. One
    GhostWorkspace is used for all samples, so nothing is allocated on the device. =#
    current = function (ctx, backend, FT)
        n = ctx.n
        domain = make_model(load_floes(0), n).domain
        src = adapt(backend, to_float(FT, Subzero.FixedWidthFloes(load_floes(n))))
        d = adapt(backend, deepcopy(src))
        workspace = Subzero.GhostWorkspace(d)
        # The ghosts must fit in the spare rows, else add_ghosts! allocates new arrays
        Subzero.add_ghosts!(d, domain, backend; workspace).poly === d.poly ||
            error("the ghosts don't fit in the capacity of the floes")
        reset!(d, src, backend)
        trial = @bench(begin
                g = Subzero.add_ghosts!($d, $domain, $backend; workspace = $workspace)
                Subzero.remove_ghosts!(g, $n, $backend)
                KernelAbstractions.synchronize($backend)
            end,
            reset!($d, $src, $backend))
        record!(ctx, "add and remove ghosts", trial; FT, primary = true)
        trial = @bench(Subzero.add_ghosts!($d, $domain, $backend; workspace = $workspace),
            reset!($d, $src, $backend))
        record!(ctx, "add_ghosts!", trial; FT)
        free!(d)
        free!(src)
        foreach(unsafe_free!, (workspace.counts, workspace.offsets))
    end,
    baseline_result = function (n)
        floes, _ = ghosted_floes(n)
        Dict("centroid" => collect(floes.centroid), "ghosts" => collect(floes.ghosts))
    end,
    current_result = function (n, backend)
        domain = make_model(load_floes(0), n).domain
        d = adapt(backend, Subzero.FixedWidthFloes(load_floes(n)))
        floes = adapt(Array, Subzero.add_ghosts!(d, domain, backend))
        rows = 1:floes.n_floes
        Dict(
            "centroid" => [floes.centroid[i, :] for i in rows],
            "ghosts" => [floes.ghosts[i, 1:floes.num_ghosts[i]] for i in rows],
        )
    end,
    float32 = true,
)
