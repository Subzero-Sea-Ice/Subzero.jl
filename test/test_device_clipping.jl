#=
The polygon clipping in `src/tools/device_clipping.jl` replaces GeometryOps' Foster-Hormann
clipping in GPU kernels, and must give the same results bit for bit. These tests compare each
step with the GeometryOps function it follows, node by node, on the CPU, and then run the
steps in a kernel on every backend.
=#
@testset "Device clipping" begin
    go_alg = GO.FosterHormannClipping(GO.Planar())
    # A GeometryOps PolyNode as the NamedTuple that _get_node returns
    _go_node(n) = (; point = n.point, inter = n.inter, neighbor = n.neighbor, idx = n.idx,
        ent_exit = n.ent_exit, crossing = n.crossing, endpoint = Int8(Int(n.endpoint)),
        fracs = n.fracs)
    _nodes(l, n) = [Subzero._get_node(l, m) for m in 1:n]
    _indices(l, n) = [l[m] for m in 1:n]

    #= Build the node lists of polygons p and q (vectors of (x, y), not closed), clipping in
    T, with GeometryOps and with the device functions. Returns whether they're the same after
    each step, and GeometryOps' lists. =#
    function _compare_build_lists(::Type{T}, p, q) where T
        ring_a, ring_b = GI.LinearRing([p; [p[1]]]), GI.LinearRing([q; [q[1]]])
        a_list, a_idx_list, n_b_intrs = GO._build_a_list(go_alg, T, ring_a, ring_b;
            exact = GO.True())
        a_after_a, a_idx_after_a = _go_node.(a_list), copy(a_idx_list)
        b_list = GO._build_b_list(go_alg, T, a_idx_list, a_list, n_b_intrs, ring_b)

        _, (poly_a, poly_b) = _array_polys([p, q])
        capacity = 2 * (length(p) + length(q)) + 64
        a = Subzero._PolyNodes(Subzero._PolyNodes(T, 1, capacity), 1)
        b = Subzero._PolyNodes(Subzero._PolyNodes(T, 1, capacity), 1)
        a_idx = Subzero._IndexList(Subzero._IndexList(1, capacity), 1)
        na, n_idx, nbi = Subzero._build_a_list!(a, a_idx, T, poly_a, poly_b)
        same_a = isequal((na, n_idx, nbi), (length(a_list), length(a_idx_after_a), n_b_intrs)) &&
            isequal(_nodes(a, na), a_after_a) && isequal(_indices(a_idx, n_idx), a_idx_after_a)
        nb = Subzero._build_b_list!(b, a, a_idx, n_idx, T, poly_b)
        same_b = nb == length(b_list) && isequal(_nodes(b, nb), _go_node.(b_list)) &&
            isequal(_nodes(a, na), _go_node.(a_list)) &&
            isequal(_indices(a_idx, n_idx), a_idx_list)
        return same_a, same_b, a_list, b_list
    end

    # Pairs that hit the special cases of building the lists, by name
    square(x0, y0, s) = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s)]
    special_pairs = (
        crossing = (square(0, 0, 2), square(1, 1, 2)),
        identical = (square(0, 0, 2), square(0, 0, 2)),
        identical_reversed = (square(0, 0, 2), reverse(square(0, 0, 2))),
        shared_edge = (square(0, 0, 2), square(2, 0, 2)),
        partly_shared_edge = (square(0, 0, 2), square(1, -2, 2)),
        touching_corners = (square(0, 0, 2), square(2, 2, 2)),
        vertex_on_edge = (square(0, 0, 2), [(1, 2), (2, 3), (0, 3)]),
        nested = (square(0, 0, 4), square(1, 1, 2)),
        nested_sharing_edge = (square(0, 0, 4), square(0, 1, 2)),
        disjoint = (square(0, 0, 1), square(3, 3, 1)),
        # collinear points along a shared edge: an overlapping chain of edges
        collinear_chain = ([(0, 0), (1, 0), (2, 0), (3, 0), (3, 2), (0, 2)],
            [(0.5, -1), (2.5, -1), (2.5, 0), (2, 0), (1, 0), (0.5, 0)]),
        # a repeated point: an edge of length zero
        repeated_point = ([(0, 0), (2, 0), (2, 0), (2, 2), (0, 2)], square(1, 1, 2)),
    )

    @testset "Building the lists: special cases with floes $FT, clipping in $T" for
            FT in (Float64, Float32), T in unique((Float64, FT))
        for (name, (p, q)) in pairs(special_pairs)
            p, q = [FT.(v) for v in p], [FT.(v) for v in q]
            for (r, s) in ((p, q), (q, p))
                same_a, same_b = _compare_build_lists(T, r, s)
                @test same_a
                @test same_b
            end
        end
    end

    #= Many pairs: the synthetic polygons with each other (including with themselves), and
    real floe shapes moved onto each other, at large coordinates too. =#
    @testset "Building the lists: random and real pairs with floes $FT, clipping in $T, \
            offset $offset" for FT in (Float64, Float32), T in unique((Float64, FT)),
            offset in (0.0, 1e6)
        rng = Xoshiro(11)
        synthetic = _geom_test_polys(Float64; n = 60)
        real_shapes = _floe_shape_rings()[1:8:end]
        test_pairs = [(synthetic[i], synthetic[j]) for i in eachindex(synthetic)
            for j in i:length(synthetic) if i == j || rand(rng) < 0.2]
        for p in real_shapes
            q = rand(rng, real_shapes)
            cp = (mean(first, p), mean(last, p))
            cq = (mean(first, q), mean(last, q))
            Δx, Δy = cp[1] - cq[1] + 300 * randn(rng), cp[2] - cq[2]
            push!(test_pairs, (p, [(x + Δx, y + Δy) for (x, y) in q]))
        end
        shift(p) = [(FT(x + offset), FT(y + offset)) for (x, y) in p]
        same = trues(2)
        n_crossing, n_touching, n_b_vertex = 0, 0, 0
        for (p, q) in test_pairs
            same_a, same_b, a_list, b_list = _compare_build_lists(T, shift(p), shift(q))
            same .&= (same_a, same_b)
            n_crossing += any(n -> n.crossing, a_list)
            # hinges and overlaps: intersection nodes that aren't crossings yet
            n_touching += any(n -> n.inter && !n.crossing, a_list)
            # B's vertices that are intersection points
            n_b_vertex += any(n -> n.inter && n.fracs[2] == 0, b_list)
        end
        @test all(same)
        @test n_crossing > 100
        @test n_touching > 50
        @test n_b_vertex > 50
    end

    @testset "Building the lists: too little room" begin
        _, (poly_a, poly_b) = _array_polys([square(0.0, 0.0, 2.0), square(1.0, 1.0, 2.0)])
        lists(capacity) = Subzero._PolyNodes(Subzero._PolyNodes(Float64, 1, capacity), 1)
        indices(capacity) = Subzero._IndexList(Subzero._IndexList(1, capacity), 1)
        # 4 vertices and 2 intersection points in each list
        @test Subzero._build_a_list!(lists(6), indices(2), Float64, poly_a, poly_b) == (6, 2, 2)
        @test Subzero._build_a_list!(lists(5), indices(2), Float64, poly_a, poly_b) == (-1, 0, 0)
        @test Subzero._build_a_list!(lists(6), indices(1), Float64, poly_a, poly_b) == (-1, 0, 0)
        a, a_idx = lists(6), indices(2)
        Subzero._build_a_list!(a, a_idx, Float64, poly_a, poly_b)
        @test Subzero._build_b_list!(lists(6), a, a_idx, 2, Float64, poly_b) == 6
        Subzero._build_a_list!(a, a_idx, Float64, poly_a, poly_b)
        @test Subzero._build_b_list!(lists(5), a, a_idx, 2, Float64, poly_b) == -1
    end

    #= Kernel test: each floe builds the lists of itself and a partner floe in its rows of
    the lists, clipping in Float64 as intersect_polys does =#
    function _build_lists_kernel!(floes, i, a_nodes, b_nodes, a_indices, counts, partner)
        a = Subzero._PolyNodes(a_nodes, i)
        b = Subzero._PolyNodes(b_nodes, i)
        a_idx = Subzero._IndexList(a_indices, i)
        poly_a, poly_b = Subzero._floe_poly(floes, i), Subzero._floe_poly(floes, partner[i])
        na, n_idx, nbi = Subzero._build_a_list!(a, a_idx, Float64, poly_a, poly_b)
        counts[i, 1], counts[i, 2], counts[i, 3] = na, n_idx, nbi
        counts[i, 4] = na < 0 ? 0 : Subzero._build_b_list!(b, a, a_idx, n_idx, Float64, poly_b)
        return
    end
    @testset "Building the lists on $backend with floes $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        polys = _geom_test_polys(FT; n = 100)
        fs = FloeSettings(FT; min_floe_area = 0)
        rng = Xoshiro(4)
        floes = StructArray([
            Floe{FT}([[[x, y] for (x, y) in [p; [p[1]]]]], FT(0.25); floe_settings = fs, rng)
            for p in polys
        ])
        fwf = Subzero.FixedWidthFloes(floes)
        n = length(floes)
        # partners include the floe itself
        partner = Int32[rand(rng) < 0.1 ? i : rand(rng, 1:n) for i in 1:n]
        capacity = 64
        outputs() = (Subzero._PolyNodes(Float64, n, capacity),
            Subzero._PolyNodes(Float64, n, capacity), Subzero._IndexList(n, capacity),
            zeros(Int32, n, 4))
        ref = outputs()
        for i in 1:n
            _build_lists_kernel!(fwf, i, ref..., partner)
        end
        dev = Subzero.adapt(backend, outputs())
        Subzero.launch_per_floe!(_build_lists_kernel!, backend, Subzero.adapt(backend, fwf),
            dev..., Subzero.adapt(backend, partner))
        Subzero.KernelAbstractions.synchronize(backend)
        got = Subzero.adapt(Array, dev)
        @test all(>=(0), ref[4])  # the capacity is enough
        @test isequal(got[4], ref[4])
        same_arrays(g, r) = all(fieldnames(typeof(r))) do name
            name == :k || isequal(getfield(g, name), getfield(r, name))
        end
        @test all(same_arrays.(got[1:3], ref[1:3]))
        # most pairs intersect
        @test count(>(0), ref[4][:, 2]) > n ÷ 2
    end
end
