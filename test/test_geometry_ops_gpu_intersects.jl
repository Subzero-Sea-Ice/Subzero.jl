#=
Tests of GeometryOpsGPU's `intersects` and `disjoint` for polygons, and the functions they use,
against GeometryOps. With Float64 polygons the results must be the same, bit for bit. Like
`test_geometry_ops_gpu.jl`, the tests use only GeometryOpsGPU, GeometryOps and the helpers in
`geometry_ops_gpu_utils.jl`.
=#
@testset "GeometryOpsGPU intersects" begin
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

end
