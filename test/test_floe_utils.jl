@testset "Floe Utils" begin
    ext = [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    hole1 = [[0.2, 0.3], [0.2, 0.2], [0.3, 0.2], [0.3, 0.3], [0.2, 0.3]]
    hole2 = [[0.5, 0.6], [0.5, 0.5], [0.6, 0.5], [0.6, 0.6], [0.5, 0.6]]
    poly_nohole = Subzero.make_polygon([ext])
    poly_hole1 = Subzero.make_polygon([ext, hole1])
    poly_hole2 = Subzero.make_polygon([ext, hole1, hole2])
    
    # Test validating/correcting RingVecs and PolyVecs
    @test Subzero.valid_ringvec!(ext) == ext
    invalid_ext = [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
    @test Subzero.valid_ringvec!(invalid_ext) == ext
    @test_throws AssertionError Subzero.valid_ringvec!([[0.0, 1.0], [0.0, 0.0]])
    invalid_coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
                      [[0.2, 0.3], [0.2, 0.2], [0.3, 0.2], [0.3, 0.3]]]
    @test Subzero.valid_polyvec!(invalid_coords) == [ext, hole1]
    duplicate_invalid_ext = [
        [0.0, 1.0],
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 0.0], 
        [1.0, 1.0],
    ]
    @test Subzero.valid_ringvec!(duplicate_invalid_ext) == ext
    @test Subzero.valid_ringvec!(
        [[0.0, 1.0], [0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
    ) == ext
    @test Subzero.valid_ringvec!(
        [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
    ) == ext
    @test_throws AssertionError Subzero.valid_ringvec!(
        [[0.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
    )
    @test_throws AssertionError Subzero.valid_polyvec!([[Float64[]]])

    # Test predicate hashole for polygons and multipolygons
    @test !Subzero.hashole(poly_nohole)
    @test Subzero.hashole(poly_hole1)
    @test Subzero.hashole(poly_hole2)

    # Test removing holes from polygons
    poly_copy_holes = Subzero.make_polygon([ext, hole1])
    Subzero.rmholes!(poly_copy_holes)
    @test GI.nhole(poly_copy_holes) == 0

    # Test moment of intertia calculations - compared to values output my MATLAB
    poly_moment = Subzero._calc_moment_inertia(Float64, Subzero.make_polygon([ext]), [0.5, 0.5], 0.25)
    @test isapprox(poly_moment, 38.333, atol = 0.001)
    @test Subzero._calc_moment_inertia(Float64, poly_nohole, GO.centroid(poly_nohole), 0.25) == poly_moment
    tri_poly = Subzero.make_polygon([[[0, 1], [0, 0], [1, 0], [0, 1]]] .* 6.67)
    tri_moment = Subzero._calc_moment_inertia(Float64, tri_poly, GO.centroid(tri_poly), 0.5)
    @test isapprox(tri_moment, 50581.145, atol = 0.001)

    @testset "_move_floe! matches _move_poly" begin
        rng = Xoshiro(1)
        for FT in (Float64, Float32)
            floes = _make_timestep_test_floes(FT, FloeSettings(FT))
            fwf = Subzero.FixedWidthFloes(floes)
            for i in eachindex(floes)
                Δx, Δy = 100 .* (rand(rng, FT, 2) .- FT(0.5))
                Δα = FT(π) * (rand(rng, FT) - FT(0.5))
                cx, cy = floes.centroid[i]
                expected = Subzero._move_poly(FT, floes.poly[i], Δx, Δy, Δα, cx, cy)
                Subzero._move_floe!(fwf, i, Δx, Δy, Δα)
                n = fwf.n_points[i]
                moved = [(fwf.poly[i, j, 1], fwf.poly[i, j, 2]) for j in 1:n]
                @test isapprox(
                    reinterpret(FT, moved),
                    reinterpret(FT, collect(GI.getpoint(expected)));
                    rtol = 10eps(FT),
                )
                @test fwf.centroid[i, :] ≈ [cx + Δx, cy + Δy]
                @test all(fwf.poly[i, (n + 1):end, :] .== 0)
            end
        end
    end
    @testset "_move_poly! matches _move_poly around the origin" begin
        for FT in (Float64, Float32)
            floes = _make_timestep_test_floes(FT, FloeSettings(FT))
            fwf = Subzero.FixedWidthFloes(floes)
            Δx, Δy, Δα = FT(-31.4), FT(27.2), FT(0.7)
            for i in eachindex(floes)
                expected = Subzero._move_poly(FT, floes.poly[i], Δx, Δy, Δα)
                Subzero._move_poly!(fwf, i, Δx, Δy, Δα)
                moved = [(fwf.poly[i, j, 1], fwf.poly[i, j, 2]) for j in 1:fwf.n_points[i]]
                @test isapprox(
                    reinterpret(FT, moved),
                    reinterpret(FT, collect(GI.getpoint(expected)));
                    rtol = 10eps(FT),
                )
            end
        end
    end
    @testset "_move_floe! on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        floes = _make_timestep_test_floes(FT, FloeSettings(FT))
        Δx, Δy, Δα = FT(-31.4), FT(27.2), FT(0.7)
        fwf = _launch_on(Subzero._move_floe!, backend, Subzero.FixedWidthFloes(floes), Δx, Δy, Δα)
        for i in eachindex(floes)
            cx, cy = floes.centroid[i]
            expected = Subzero._move_poly(FT, floes.poly[i], Δx, Δy, Δα, cx, cy)
            n = fwf.n_points[i]
            moved = [(fwf.poly[i, j, 1], fwf.poly[i, j, 2]) for j in 1:n]
            @test isapprox(
                reinterpret(FT, moved),
                reinterpret(FT, collect(GI.getpoint(expected)));
                rtol = 10eps(FT),
            )
            @test fwf.centroid[i, :] ≈ [cx + Δx, cy + Δy]
            @test all(fwf.poly[i, (n + 1):end, :] .== 0)
        end
    end
    #=
    Floes for testing `_overlaps_rectangle`, and rectangles to test them with. Besides
    random star-shaped floes, there are axis-aligned rectangles and L-shapes on a grid with a
    spacing of 250 m, and the rectangles' sides are on the same grid or go through floe
    vertices, so floes and rectangles often share edges or vertices. The floes from
    `edge_cases` are tested with the rectangle 10 ≤ x ≤ 15 km, -5 ≤ y ≤ 15 km, the eastern
    boundary of a 10 km × 10 km domain.
    =#
    function _overlap_test_floes(::Type{FT}; n = 200) where FT
        rng = Xoshiro(3)
        fs = FloeSettings(FT; min_floe_area = 0)
        edge_cases = [
            [(1, 1), (3, 1), (3, 3), (1, 3)],  # inside the domain
            [(8, 1), (10, 1), (10, 3), (8, 3)],  # edge on the rectangle's side
            [(8, 2), (10, 3), (8, 4)],  # vertex on the rectangle's side
            [(8, 1), (11, 1), (11, 3), (8, 3)],  # crosses the side
            [(8, 1), (11, 1), (11, 3), (8, 3)][end:-1:1],  # the same, clockwise
            [(8, 2), (10.5, 3), (8, 4)],  # one vertex inside
            [(11, 1), (12, 1), (12, 3), (11, 3)],  # entirely inside
            [(15, 1), (17, 1), (17, 3), (15, 3)],  # edge on the far side
            [(16, 1), (17, 1), (17, 3), (16, 3)],  # beyond the far side
            [(5, -10), (20, -10), (20, 20), (5, 20)],  # contains the rectangle
            [(9, -6), (20, -6), (9, -4)],  # an edge crosses, no vertex inside
            [(8, -7), (10, -5), (8, -6)],  # vertex on the corner
        ]
        polys = [[(FT(1e3x), FT(1e3y)) for (x, y) in pts] for pts in edge_cases]
        for _ in 1:n
            if rand(rng) < 0.5
                nv = rand(rng, 3:9)
                θs = sort(2π .* rand(rng, nv))
                rs = rand(rng) < 0.5 ? fill(1e3, nv) : 1e3 .* (0.3 .+ rand(rng, nv))
                push!(polys, [(FT(r * cos(θ)), FT(r * sin(θ))) for (r, θ) in zip(rs, θs)])
            else
                g = FT.(250 .* sort(rand(rng, -4:4, 3)))
                h = FT.(250 .* sort(rand(rng, -4:4, 3)))
                (allunique(g) && allunique(h)) || continue
                push!(polys, rand(rng) < 0.5 ?
                    [(g[1], h[1]), (g[3], h[1]), (g[3], h[3]), (g[1], h[3])] :
                    [(g[1], h[1]), (g[3], h[1]), (g[3], h[2]), (g[2], h[2]), (g[2], h[3]), (g[1], h[3])])
            end
        end
        floes = StructArray([
            Floe{FT}([[[x, y] for (x, y) in [p; [p[1]]]]], FT(0.25); floe_settings = fs, rng)
            for p in polys
        ])
        vertices = [p for poly in polys[(length(edge_cases) + 1):end] for p in poly]
        pick(i) = (u = rand(rng); u < 0.4 ? rand(rng, vertices)[i] :
            u < 0.7 ? FT(250 * rand(rng, -4:4)) : FT(4e3 * rand(rng) - 2e3))
        rects = NTuple{4, FT}[]
        while length(rects) < 20
            xmin, xmax = minmax(pick(1), pick(1))
            ymin, ymax = minmax(pick(2), pick(2))
            xmin < xmax && ymin < ymax && push!(rects, (xmin, xmax, ymin, ymax))
        end
        return floes, length(edge_cases), rects
    end
    overlaps_go(floes, i, (xmin, xmax, ymin, ymax), FT) = !isempty(Subzero.intersect_polys(
        floes.poly[i], Subzero._make_bounding_box_polygon(FT, xmin, xmax, ymin, ymax), FT))
    @testset "_overlaps_rectangle matches intersect_polys with $FT" for FT in (Float64, Float32)
        floes, n_edge_cases, rects = _overlap_test_floes(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        bound = PeriodicBoundary(East, FT; x0 = 0, xf = 1e4, y0 = 0, yf = 1e4)
        ext = GI.extent(bound.poly)
        east = (ext.X..., ext.Y...)
        @test east == (FT(1e4), FT(1.5e4), FT(-5e3), FT(1.5e4))
        @test [Subzero._overlaps_rectangle(fwf, i, east...) for i in 1:n_edge_cases] ==
            [false, false, false, true, true, true, true, false, false, true, true, false]
        for rect in [east; rects], i in eachindex(floes)
            @test Subzero._overlaps_rectangle(fwf, i, rect...) == overlaps_go(floes, i, rect, FT)
        end
    end
    #= Kernel test: mark_overlap! writes whether each floe overlaps the rectangle into its
    flags. These are test markers, not FLAG_* event bits. =#
    mark_overlap!(floes, i, rect...) =
        (floes.flags[i] = Subzero._overlaps_rectangle(floes, i, rect...) ? 0x01 : 0x00; nothing)
    @testset "_overlaps_rectangle on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        floes, _, rects = _overlap_test_floes(FT)
        for rect in rects
            fwf = _launch_on(mark_overlap!, backend, Subzero.FixedWidthFloes(floes), rect...)
            @test fwf.flags[1:length(floes)] ==
                [overlaps_go(floes, i, rect, FT) ? 0x01 : 0x00 for i in eachindex(floes)]
        end
    end
end