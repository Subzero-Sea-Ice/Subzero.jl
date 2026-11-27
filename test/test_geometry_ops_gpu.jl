#=
Tests of the GeometryOpsGPU submodule against GeometryOps (through Subzero's CPU geometry
functions, which call it). With Float64 polygons the results must be the same, bit for bit,
except for `RotateTranslate`, which matches to round-off; with Float32 polygons they may
differ by round-off.
=#
@testset "GeometryOpsGPU" begin
    GOG = Subzero.GeometryOpsGPU
    Extent = Subzero.Extents.Extent
    # The points of a GeoInterface polygon's exterior ring, as tuples
    gi_points(p) = [(GI.x(q), GI.y(q)) for q in GI.getpoint(GI.getexterior(p))]

    @testset "ArrayPolygon" begin
        coords = zeros(2, 4, 2)
        coords[2, 1:4, 1] .= [0, 1, 1, 0]
        coords[2, 1:4, 2] .= [0, 0, 1, 0]
        p = GOG.ArrayPolygon(coords, 2, 4)
        @test GOG.npoint(p) == 4
        @test [GOG.getpoint(p, j) for j in 1:4] == [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]
        # row 1 is left alone
        GOG.transform!(GOG.Translate(1.0, 2.0), p)
        @test [GOG.getpoint(p, j) for j in 1:4] == [(1.0, 2.0), (2.0, 2.0), (2.0, 3.0), (1.0, 2.0)]
        @test all(iszero, coords[1, :, :])
        # it's a GeoInterface polygon without holes, which GeometryOps accepts
        @test GI.testgeometry(p)
        @test GI.nhole(p) == 0
        @test gi_points(p) == [GOG.getpoint(p, j) for j in 1:4]
        @test GO.area(p) == 0.5
    end

    # transform! changes the polygon in place and leaves the padding after its points alone
    @testset "transform! with $FT" for FT in (Float64, Float32)
        floes = _make_timestep_test_floes(FT, FloeSettings(FT))
        n = length(floes)
        # Translate adds Δx and Δy to every coordinate, as GO.transform with a Translation
        Δx, Δy = FT(-2e4), FT(1234.5)
        fwf = Subzero.FixedWidthFloes(floes)
        foreach(i -> GOG.transform!(GOG.Translate(Δx, Δy), Subzero._floe_poly(fwf, i)), 1:n)
        expected = [Subzero._translate_poly(FT, floes.poly[i], Δx, Δy) for i in 1:n]
        @test isempty(_failing(i -> _fwf_points(fwf, i) == collect(GI.getpoint(expected[i])), 1:n))
        @test isempty(_failing(i -> _fwf_padding_is_zero(fwf, i), 1:n))
        # RotateTranslate, as GO.transform with Translation ∘ recenter(LinearMap(Angle2d)),
        # which computes in another order, so up to round-off
        Δα = FT(0.7)
        fwf = Subzero.FixedWidthFloes(floes)
        foreach(1:n) do i
            cx, cy = floes.centroid[i]
            GOG.transform!(GOG.RotateTranslate(Δα, cx, cy, Δx, Δy), Subzero._floe_poly(fwf, i))
        end
        expected = [Subzero._move_poly(FT, floes.poly[i], Δx, Δy, Δα, floes.centroid[i]...)
            for i in 1:n]
        @test isempty(_not_approx(
            [reinterpret(FT, _fwf_points(fwf, i)) for i in 1:n],
            [reinterpret(FT, collect(GI.getpoint(p))) for p in expected];
            rtol = 10eps(FT),
        ))
        @test isempty(_failing(i -> _fwf_padding_is_zero(fwf, i), 1:n))
    end

    #=
    Floes for testing `interiors_intersect`, and rectangles to test them with. Besides
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
        rects = Extent[]
        while length(rects) < 20
            xmin, xmax = minmax(pick(1), pick(1))
            ymin, ymax = minmax(pick(2), pick(2))
            xmin < xmax && ymin < ymax && push!(rects, Extent(X = (xmin, xmax), Y = (ymin, ymax)))
        end
        return floes, length(edge_cases), rects
    end
    overlaps_go(floes, i, ext, FT) = !isempty(Subzero.intersect_polys(floes.poly[i],
        Subzero._make_bounding_box_polygon(FT, ext.X..., ext.Y...), FT))
    intersects(fwf, i, ext) = GOG.interiors_intersect(Subzero._floe_poly(fwf, i), ext)
    @testset "interiors_intersect matches intersect_polys with $FT" for FT in (Float64, Float32)
        floes, n_edge_cases, rects = _overlap_test_floes(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        east = GI.extent(PeriodicBoundary(East, FT; x0 = 0, xf = 1e4, y0 = 0, yf = 1e4).poly)
        @test (east.X..., east.Y...) == (FT(1e4), FT(1.5e4), FT(-5e3), FT(1.5e4))
        @test [intersects(fwf, i, east) for i in 1:n_edge_cases] ==
            [false, false, false, true, true, true, true, false, false, true, true, false]
        # (rectangle, floe) pairs where the two differ. In Float32 they may differ by
        # round-off, but don't for these floes.
        mismatches = [(rect, i) for rect in [east; rects], i in eachindex(floes)
            if intersects(fwf, i, rect) != overlaps_go(floes, i, rect, FT)]
        @test isempty(mismatches)
    end
    # Bounds of another type than the polygon's aren't accepted
    @testset "interiors_intersect with Float32 polygons" begin
        floes, _, rects = _overlap_test_floes(Float32)
        fwf = Subzero.FixedWidthFloes(floes)
        @test_throws MethodError intersects(fwf, 1, Extent(X = (1e4, 1.5e4), Y = (-5e3, 1.5e4)))
    end
    #= Kernel test: mark_overlap! writes whether each floe overlaps the rectangle into its
    flags. These are test markers, not FLAG_* event bits. It names the module in full: the
    local `GOG` would be captured, and a module can't be passed to a kernel. =#
    mark_overlap!(floes, i, ext) = (floes.flags[i] =
        Subzero.GeometryOpsGPU.interiors_intersect(Subzero._floe_poly(floes, i), ext) ? 0x01 : 0x00;
        nothing)
    @testset "interiors_intersect on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        floes, _, rects = _overlap_test_floes(FT)
        marks(ext) = _launch_on(mark_overlap!, backend, Subzero.FixedWidthFloes(floes), ext).flags[1:length(floes)]
        expected(ext) = [overlaps_go(floes, i, ext, FT) ? 0x01 : 0x00 for i in eachindex(floes)]
        # rectangles for which the kernel's result differs
        @test isempty(_failing(ext -> marks(ext) == expected(ext), rects))
    end
end
