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
end
