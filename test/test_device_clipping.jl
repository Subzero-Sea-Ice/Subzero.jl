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

    # The regions in `out` and the polygons that GeometryOps returns as vectors of points
    _regions(out, n) = [[Subzero.GeometryOpsGPU.getpoint(Subzero._region(out, m), j)
        for j in 1:Subzero.GeometryOpsGPU.npoint(Subzero._region(out, m))] for m in 1:n]
    _go_regions(polys) = [collect(GI.getpoint(GI.getexterior(r))) for r in polys]

    #= Clip polygons p and q (vectors of (x, y), not closed) in T with GeometryOps and with
    the device functions. Returns whether the node lists after all steps of
    `GO._build_ab_list` are the same, whether the traced regions are the same, whether the
    result of `intersect_polys!` is the same as `intersect_polys` (only for T = Float64, the
    type that `intersect_polys` clips in), and what the pair tests (for counting cases). =#
    function _compare_clipping(::Type{T}, p, q) where T
        ring_a, ring_b = GI.LinearRing([p; [p[1]]]), GI.LinearRing([q; [q[1]]])
        a_list, b_list, a_idx_list = GO._build_ab_list(go_alg, T, ring_a, ring_b,
            GO._inter_delay_cross_f, GO._inter_delay_bounce_f; exact = GO.True())
        traced = GO._trace_polynodes(go_alg, T, a_list, b_list, copy(a_idx_list),
            GO._inter_step, ring_a, ring_b)

        _, (poly_a, poly_b) = _array_polys([p, q])
        capacity = 2 * (length(p) + length(q)) + 64
        lists() = Subzero._ClipLists(Subzero._ClipLists(T, 1, capacity), 1)
        regions() = Subzero._PolyRegions(Subzero._PolyRegions(T, 1, 16, capacity), 1)
        l = lists()
        na, n_idx, _ = Subzero._build_a_list!(l.a, l.a_idx, T, poly_a, poly_b)
        nb = Subzero._build_b_list!(l.b, l.a, l.a_idx, n_idx, T, poly_b)
        ok = Subzero._classify_crossing!(l.a, na, l.b, nb) &&
            Subzero._flag_ent_exit!(l.a, na, poly_b) && Subzero._flag_ent_exit!(l.b, nb, poly_a)
        n_cross = Subzero._index_crossing_intrs!(l.a, l.b, l.a_idx, n_idx)
        same_lists = ok && isequal(_nodes(l.a, na), _go_node.(a_list)) &&
            isequal(_nodes(l.b, nb), _go_node.(b_list)) &&
            isequal(_indices(l.a_idx, n_cross), a_idx_list)
        out = regions()
        n_traced = Subzero._trace_polynodes!(out, l.a, na, l.b, nb, l.a_idx, n_cross)
        same_traced = isequal(_regions(out, n_traced), _go_regions(traced))

        result = nothing
        same_result = true
        if T == Float64
            FT = eltype(eltype(p))
            result = _go_regions(Subzero.intersect_polys(_gi_poly(p), _gi_poly(q), FT))
            out = regions()
            n = Subzero.intersect_polys!(out, lists(), poly_a, poly_b)
            same_result = n >= 0 && isequal(_regions(out, n), result)
        end
        cases = (;
            chain = any(n -> n.endpoint != GO.not_endpoint, a_list),
            # the boundaries don't cross, but one polygon is inside the other
            nested = isempty(traced) && !isnothing(result) && !isempty(result),
            several = !isnothing(result) && length(result) > 1,
            collinear_removed = !isnothing(result) &&
                sum(length, result; init = 0) < sum(r -> GI.npoint(r), traced; init = 0),
        )
        return same_lists, same_traced, same_result, cases
    end

    @testset "Clipping: special cases with floes $FT, clipping in $T" for
            FT in (Float64, Float32), T in unique((Float64, FT))
        for (name, (p, q)) in pairs(special_pairs)
            p, q = [FT.(v) for v in p], [FT.(v) for v in q]
            for (r, s) in ((p, q), (q, p))
                same_lists, same_traced, same_result = _compare_clipping(T, r, s)
                @test same_lists
                @test same_traced
                @test same_result
            end
        end
    end

    #= The same pairs as for building the lists, but every pair of synthetic polygons, so
    that all the special cases come up many times =#
    @testset "Clipping: random and real pairs with floes $FT, clipping in $T, \
            offset $offset" for FT in (Float64, Float32), T in unique((Float64, FT)),
            offset in (0.0, 1e6)
        rng = Xoshiro(11)
        synthetic = _geom_test_polys(Float64; n = 60)
        real_shapes = _floe_shape_rings()[1:8:end]
        test_pairs = [(synthetic[i], synthetic[j]) for i in eachindex(synthetic)
            for j in i:length(synthetic)]
        for p in real_shapes
            q = rand(rng, real_shapes)
            cp = (mean(first, p), mean(last, p))
            cq = (mean(first, q), mean(last, q))
            Δx, Δy = cp[1] - cq[1] + 300 * randn(rng), cp[2] - cq[2]
            push!(test_pairs, (p, [(x + Δx, y + Δy) for (x, y) in q]))
        end
        shift(p) = [(FT(x + offset), FT(y + offset)) for (x, y) in p]
        same = trues(3)
        counts = Dict(:chain => 0, :nested => 0, :several => 0, :collinear_removed => 0)
        for (p, q) in test_pairs, (r, s) in ((shift(p), shift(q)), (shift(q), shift(p)))
            same_lists, same_traced, same_result, cases = _compare_clipping(T, r, s)
            same .&= (same_lists, same_traced, same_result)
            for (case, happened) in pairs(cases)
                counts[case] += happened
            end
        end
        @test all(same)
        @test counts[:chain] > 100
        if T == Float64
            @test counts[:nested] > 100
            @test counts[:several] > 50
            @test counts[:collinear_removed] > 30
        end
    end

    @testset "Clipping: too little room" begin
        polys = [square(0.0, 0.0, 2.0), square(1.0, 1.0, 2.0), square(0.5, 0.5, 1.0),
            # a U shape and a bar across both of its arms: two regions
            [(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (2.0, 3.0), (2.0, 1.0), (1.0, 1.0),
                (1.0, 3.0), (0.0, 3.0)],
            [(-1.0, 2.0), (4.0, 2.0), (4.0, 2.5), (-1.0, 2.5)]]
        _, (crossing_a, crossing_b, inside, u_shape, bar) = _array_polys(polys)
        lists(capacity) = Subzero._ClipLists(Subzero._ClipLists(Float64, 1, capacity), 1)
        regions(max_polys, capacity) =
            Subzero._PolyRegions(Subzero._PolyRegions(Float64, 1, max_polys, capacity), 1)
        clip(l, out, p, q) = Subzero.intersect_polys!(out, l, p, q)
        too_small = Subzero._CLIP_TOO_SMALL
        # crossing squares: 6 nodes in each list, and a region of 4 points and the
        # repeated first point
        @test clip(lists(6), regions(1, 5), crossing_a, crossing_b) == 1
        @test clip(lists(5), regions(1, 5), crossing_a, crossing_b) == too_small
        @test clip(lists(6), regions(1, 4), crossing_a, crossing_b) == too_small
        @test clip(lists(6), regions(0, 5), crossing_a, crossing_b) == too_small
        # a square inside the other: the region is a copy of it
        @test clip(lists(8), regions(1, 5), crossing_a, inside) == 1
        @test clip(lists(8), regions(1, 4), crossing_a, inside) == too_small
        @test clip(lists(8), regions(0, 5), crossing_a, inside) == too_small
        @test clip(lists(32), regions(2, 16), u_shape, bar) == 2
        @test clip(lists(32), regions(1, 16), u_shape, bar) == too_small
    end

    #= Kernel test: each floe intersects itself with a partner floe, with lists and regions
    in its rows, as the collisions and other steps will =#
    function _intersect_kernel!(floes, i, lists, regions, counts, partner)
        poly_a, poly_b = Subzero._floe_poly(floes, i), Subzero._floe_poly(floes, partner[i])
        counts[i] = Subzero.intersect_polys!(Subzero._PolyRegions(regions, i),
            Subzero._ClipLists(lists, i), poly_a, poly_b)
        return
    end
    @testset "Clipping on $backend with floes $FT" for backend in test_backends(),
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
        capacity, max_polys = 64, 8
        outputs() = (Subzero._ClipLists(Float64, n, capacity),
            Subzero._PolyRegions(Float64, n, max_polys, capacity), zeros(Int32, n))
        ref = outputs()
        for i in 1:n
            _intersect_kernel!(fwf, i, ref..., partner)
        end
        # the reference is the same as intersect_polys on the floes' polygons
        @test all(>=(0), ref[3])  # the capacity is enough
        @test all(1:n) do i
            expected = _go_regions(Subzero.intersect_polys(floes.poly[i],
                floes.poly[partner[i]], FT))
            isequal(_regions(Subzero._PolyRegions(ref[2], i), ref[3][i]), expected)
        end
        dev = Subzero.adapt(backend, outputs())
        Subzero.launch_per_floe!(_intersect_kernel!, backend, Subzero.adapt(backend, fwf),
            dev..., Subzero.adapt(backend, partner))
        Subzero.KernelAbstractions.synchronize(backend)
        got = Subzero.adapt(Array, dev)
        @test isequal(got[3], ref[3])
        @test all(1:n) do i
            isequal(_regions(Subzero._PolyRegions(got[2], i), got[3][i]),
                _regions(Subzero._PolyRegions(ref[2], i), ref[3][i]))
        end
        # most pairs overlap, and some in several regions
        @test count(>(0), ref[3]) > n ÷ 2
        @test count(>(1), ref[3]) > 0
    end
end
