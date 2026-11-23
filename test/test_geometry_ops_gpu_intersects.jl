#=
Tests of GeometryOpsGPU's `intersects` and `disjoint` for polygons, and the functions they use,
against GeometryOps. With Float64 polygons the results must be the same, bit for bit. Like
`test_geometry_ops_gpu.jl`, the tests use only GeometryOpsGPU, GeometryOps and the helpers in
`geometry_ops_gpu_utils.jl`.
=#
@testset "GeometryOpsGPU intersects" begin
    #=
    Pairs of polygons (vectors of points, not closed) with all the ways two polygons can meet:
    named special cases, in both orders; polygons that touch exactly (shifted by their width
    along an axis-aligned side) or are a few millimetres apart; every pair of the synthetic
    polygons, including each with itself; and real floe shapes moved onto each other. All
    shifted by `offset`, where larger coordinates make the floats further apart.
    =#
    function _intersects_test_pairs(::Type{FT}; offset = 0.0) where FT
        sq(x0, y0, s) = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s)]
        special = [
            (sq(0, 0, 2), sq(1, 1, 2)),  # crossing
            (sq(0, 0, 2), sq(0, 0, 2)),  # identical
            (sq(0, 0, 2), reverse(sq(0, 0, 2))),  # identical, other direction
            (sq(0, 0, 2), sq(2, 0, 2)),  # shared edge
            (sq(0, 0, 2), sq(1, -2, 2)),  # partly shared edge
            (sq(0, 0, 2), sq(2, 2, 2)),  # touching corners
            (sq(0, 0, 2), [(1, 2), (2, 3), (0, 3)]),  # vertex on an edge
            (sq(0, 0, 4), sq(1, 1, 2)),  # nested
            (sq(0, 0, 4), sq(0, 1, 2)),  # nested, sharing an edge
            (sq(0, 0, 1), sq(3, 3, 1)),  # disjoint
            (sq(0, 0, 1), sq(0.5, 1.5, 1)),  # disjoint, overlapping in x
            # collinear points along a shared edge: an overlapping chain of edges
            ([(0, 0), (1, 0), (2, 0), (3, 0), (3, 2), (0, 2)],
                [(0.5, -1), (2.5, -1), (2.5, 0), (2, 0), (1, 0), (0.5, 0)]),
            # a repeated point: an edge of length zero
            ([(0, 0), (2, 0), (2, 0), (2, 2), (0, 2)], sq(1, 1, 2)),
        ]
        pairs = [x for (p, q) in special for x in ((p, q), (q, p))]
        rng = Xoshiro(13)
        synthetic = _geom_test_polys(Float64; n = 40)
        for p in synthetic
            xs, ys = extrema(first, p), extrema(last, p)
            # touching: shifted by its width or height; and a few millimetres further
            for (Δx, Δy) in ((xs[2] - xs[1], 0), (0, ys[2] - ys[1]))
                push!(pairs, (p, [(x + Δx, y + Δy) for (x, y) in p]))
                g = 1e-3 * rand(rng, 1:5)
                push!(pairs, (p, [(x + Δx + sign(Δx) * g, y + Δy + sign(Δy) * g) for (x, y) in p]))
            end
        end
        append!(pairs, [(synthetic[i], synthetic[j]) for i in eachindex(synthetic)
            for j in i:length(synthetic)])
        real_shapes = _floe_shape_rings()[1:8:end]
        for p in real_shapes
            q = rand(rng, real_shapes)
            cp = (mean(first, p), mean(last, p))
            cq = (mean(first, q), mean(last, q))
            Δx, Δy = cp[1] - cq[1] + 300 * randn(rng), cp[2] - cq[2]
            push!(pairs, (p, [(x + Δx, y + Δy) for (x, y) in q]))
        end
        shift(p) = [(FT(x + offset), FT(y + offset)) for (x, y) in p]
        return [(shift(p), shift(q)) for (p, q) in pairs]
    end

    # All pairs of edges of the polygons, as (a1, a2, b1, b2)
    _edge_pairs(polys) = [(p[a], p[mod1(a + 1, end)], q[b], q[mod1(b + 1, end)])
        for p in polys for q in polys for a in eachindex(p) for b in eachindex(q)]

    #= Points tested against segments: the segments' ends and midpoints, which are on them
    (exactly, for the axis-aligned edges of the grid polygons), and other test points =#
    @testset "Point and segment orientations with $FT" for FT in (Float64, Float32)
        polys = _geom_test_polys(FT; n = 60)
        rng = Xoshiro(21)
        same, orientations = true, GOG.PointOrientation[]
        for _ in 1:20000
            p = rand(rng, polys)
            k = rand(rng, eachindex(p))
            s1, s2 = p[k], p[mod1(k + 1, end)]
            u = rand(rng)
            pt = u < 0.25 ? s1 : u < 0.5 ? (s1 .+ s2) ./ 2 : _geom_test_points(rng, p, polys)
            o = GOG._point_segment_orientation(pt, s1, s2)
            push!(orientations, o)
            same &= Int(o) == Int(GO._point_segment_orientation(GO.Planar(), pt, s1, s2))
        end
        @test same
        @test all(o -> count(==(o), orientations) > 1000, (GOG.point_in, GOG.point_on, GOG.point_out))
    end

    #= GeometryOps computes the segment-segment orientation in Float64. With Float64 points all
    of it must be the same; with Float32 points, the kind of intersection must be the same
    with exact orientation tests, which decide it exactly in either type. =#
    @testset "Segment orientations with $FT" for FT in (Float64, Float32)
        edges = _edge_pairs(_geom_test_polys(FT; n = 40))
        go(e, exact) = GO._seg_seg_orientation(GO.Planar(), e...; exact = exact ? GO.True() : GO.False())
        gog(e, exact) = GOG._seg_seg_orientation(e...; exact = exact ? GOG.True() : GOG.False())
        if FT == Float64
            @test all(e -> isequal((Int(gog(e, true)[1]), gog(e, true)[2:3]...),
                (Int(go(e, true)[1]), go(e, true)[2:3]...)), edges)
            @test all(e -> isequal((Int(gog(e, false)[1]), gog(e, false)[2:3]...),
                (Int(go(e, false)[1]), go(e, false)[2:3]...)), edges)
        else
            @test all(e -> Int(gog(e, true)[1]) == Int(go(e, true)[1]), edges)
        end
    end

    #= Segments that hinge with a curve (meet it at a vertex, or along an edge), which is
    when _line_filled_curve_interactions splits them: the edges of each polygon against the
    ring of each other polygon, if some edge pair hinges. Compared with GeometryOps'
    _split_segment_interactions for each combination of incoming in_curve and out_curve.

    With Float32 points and exact orientation tests, the midpoint of a piece that lies along
    an edge of the curve is on that edge for GeometryOps, which computes it exactly in
    Float64, but GOG's Float32 midpoint is rounded off the edge, so the piece is in or out.
    That is allowed where a piece lies along the curve: its midpoint, computed exactly in
    Float64, is on the curve. (For disjoint and intersects it doesn't matter: a segment is
    split only when it meets the curve, which already decides them.) =#
    @testset "Splitting segments with $FT, exact = $exact" for FT in (Float64, Float32),
            exact in (true, false)
        polys = _geom_test_polys(FT; n = 60)
        _, apolys = _array_polys(polys)
        rings = [GI.LinearRing([p; [p[1]]]) for p in polys]
        gog_exact, go_exact = exact ? (GOG.True(), GO.True()) : (GOG.False(), GO.False())
        hinges(s1, s2, q) = any(eachindex(q)) do m
            GOG._seg_seg_orientation(s1, s2, q[m], q[mod1(m + 1, end)]; exact = gog_exact)[1] ==
                GOG.line_hinge
        end
        # GOG's pieces of the segment, as (start, end) pairs
        function pieces(s1, s2, r)
            ends, has, key = [s1], false, (zero(FT), s1)
            while true
                found, key = GOG._next_segment_intersection(s1, s2, r, has, key)
                found || break
                push!(ends, key[2])
                has = true
            end
            push!(ends, s2)
            return zip(ends[1:(end - 1)], ends[2:end])
        end
        along_curve(s1, s2, j) = any(pieces(s1, s2, apolys[j])) do (a, b)
            mid = ((Float64(a[1]) + b[1]) / 2, (Float64(a[2]) + b[2]) / 2)
            GO._point_filled_curve_orientation(GO.Planar(), mid, rings[j]; exact = GO.True()) ==
                GO.point_on
        end
        same, n, n_several = true, 0, 0
        for i in eachindex(polys), j in eachindex(polys), k in eachindex(polys[i])
            s1, s2 = polys[i][k], polys[i][mod1(k + 1, end)]
            hinges(s1, s2, polys[j]) || continue
            n += 1
            n_several += length(GO.intersection_points(GI.Line([s1, s2]), rings[j])) > 1
            for flags in ((false, false), (true, false), (false, true))
                same &= GOG._split_segment_interactions(s1, s2, apolys[j], flags...;
                    exact = gog_exact) ==
                    GO._split_segment_interactions(GO.Planar(), s1, s2, rings[j], flags...;
                        exact = go_exact) ||
                    (FT == Float32 && exact && along_curve(s1, s2, j))
            end
        end
        @test same
        @test n > 1000
        # segments that meet the curve at several points, so the order matters
        @test n_several > 500
    end


    #= Where the ring of each polygon of a pair is relative to the other polygon. With
    Float32 points, the in and out of pieces along the other ring's boundary depend on how
    their midpoints are rounded (see "Splitting segments"), so they may differ where both
    find the ring on the other's boundary. =#
    @testset "Line and curve interactions with $FT, offset $offset, exact = $exact" for
            FT in (Float64, Float32), offset in (0.0, 1e6), exact in (true, false)
        gog_exact, go_exact = exact ? (GOG.True(), GO.True()) : (GOG.False(), GO.False())
        same, kinds = true, Set{NTuple{3, Bool}}()
        for (p, q) in _intersects_test_pairs(FT; offset)
            _, (a, b) = _array_polys([p, q])
            x = GOG._line_filled_curve_interactions(a, b; exact = gog_exact, closed_line = true)
            y = GO._line_filled_curve_interactions(GO.Planar(), GI.LinearRing([p; [p[1]]]),
                GI.LinearRing([q; [q[1]]]); exact = go_exact, closed_line = true)
            push!(kinds, x)
            same &= x == y || (FT == Float32 && x[2] && y[2])
        end
        @test same
        # the rings are in, on and out of the other polygons in many combinations
        @test length(kinds) >= 5
    end


    #= intersects and disjoint against GeometryOps, on all the test pairs. Unlike the line
    and curve interactions, they match in Float32 too on these pairs: wherever rounding
    changes an in or out, the polygons touch, which decides both. =#
    @testset "intersects and disjoint with $FT, offset $offset" for FT in (Float64, Float32),
            offset in (0.0, 1e6)
        same, n_apart, n_touching, n_overlapping = true, 0, 0, 0
        for (p, q) in _intersects_test_pairs(FT; offset)
            _, (a, b) = _array_polys([p, q])
            gp, gq = _gi_poly(p), _gi_poly(q)
            x = GOG.intersects(a, b)
            same &= x == GO.intersects(gp, gq) && GOG.disjoint(a, b) == GO.disjoint(gp, gq) &&
                x == !GOG.disjoint(a, b)
            if !x
                n_apart += 1
            elseif isempty(GO.intersection(gp, gq; target = GI.PolygonTrait(),
                    fix_multipoly = nothing))
                n_touching += 1
            else
                n_overlapping += 1
            end
        end
        @test same
        # pairs that are apart, only touch, and overlap
        @test n_apart > 100
        @test n_touching > 30
        @test n_overlapping > 100
    end

    # Polygons of different float types aren't accepted
    @testset "intersects with polygons of different types" begin
        _, (a,) = _array_polys([[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]])
        _, (b,) = _array_polys([[(0f0, 0f0), (1f0, 0f0), (0f0, 1f0)]])
        @test_throws MethodError GOG.intersects(a, b)
    end

    # Kernel test: whether each pair of polygons intersects, and whether it is disjoint
    function _intersects_kernel!(i, coords_a, ns_a, coords_b, ns_b, out)
        a = GOG.ArrayPolygon(coords_a, i, ns_a[i])
        b = GOG.ArrayPolygon(coords_b, i, ns_b[i])
        out[i, 1] = GOG.intersects(a, b)
        out[i, 2] = GOG.disjoint(a, b)
        return
    end
    @testset "intersects on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        pairs = _intersects_test_pairs(FT)
        coords_a, polys_a = _array_polys(first.(pairs))
        coords_b, polys_b = _array_polys(last.(pairs))
        n = length(pairs)
        expected = hcat([GOG.intersects(polys_a[i], polys_b[i]) for i in 1:n],
            [GOG.disjoint(polys_a[i], polys_b[i]) for i in 1:n])
        dev = KA.zeros(backend, Bool, n, 2)
        _launch_float32_checked(backend, FT) do
            _gog_launch!(_intersects_kernel!, backend, n, Adapt.adapt(backend, coords_a),
                Adapt.adapt(backend, GOG.npoint.(polys_a)), Adapt.adapt(backend, coords_b),
                Adapt.adapt(backend, GOG.npoint.(polys_b)), dev)
        end
        @test Array(dev) == expected
        @test count(expected[:, 1]) > 100 && count(expected[:, 2]) > 100
    end
end
