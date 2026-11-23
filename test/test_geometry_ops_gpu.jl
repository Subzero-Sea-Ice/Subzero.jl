#=
Tests of the GeometryOpsGPU submodule against GeometryOps, which accepts `ArrayPolygon`s as
GeoInterface polygons. With Float64 polygons the results must be the same, bit for bit; with
Float32 polygons they may differ by round-off.
=#
@testset "GeometryOpsGPU" begin
    GOG = Subzero.GeometryOpsGPU
    Extent = Subzero.Extents.Extent
    CT = Subzero.CoordinateTransformations
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

    #= transform! changes the polygon in place, as GO.transform with CoordinateTransformations'
    types, and leaves the padding after its points alone =#
    @testset "transform! with $FT" for FT in (Float64, Float32)
        floes = _make_timestep_test_floes(FT, FloeSettings(FT))
        n = length(floes)
        original = Subzero.FixedWidthFloes(floes)
        Δx, Δy, Δα = FT(-2e4), FT(1234.5), FT(0.7)
        cs = [floes.centroid[i] for i in 1:n]
        transforms = (
            translate = (i -> GOG.Translate(Δx, Δy), i -> CT.Translation(Δx, Δy)),
            # the centre as a tuple, as in _move_poly: with a Vector, StaticArrays multiplies
            # with muladd, which may or may not be fused, depending on the CPU
            rotate_translate = (i -> GOG.RotateTranslate(Δα, cs[i]..., Δx, Δy),
                i -> CT.Translation(Δx, Δy) ∘
                    CT.recenter(CT.LinearMap(Subzero.Rotations.Angle2d(Δα)), Tuple(cs[i]))),
        )
        @testset "$name" for (name, (gog, ct)) in pairs(transforms)
            fwf = Subzero.FixedWidthFloes(floes)
            foreach(i -> GOG.transform!(gog(i), Subzero._floe_poly(fwf, i)), 1:n)
            expected = [gi_points(GO.transform(ct(i), Subzero._floe_poly(original, i))) for i in 1:n]
            # bit for bit in Float64, to round-off in Float32
            rtol = FT == Float64 ? 0.0 : 10eps(FT)
            @test isempty(_not_approx([reinterpret(FT, _fwf_points(fwf, i)) for i in 1:n],
                [reinterpret(FT, FT.(collect(Iterators.flatten(e)))) for e in expected]; rtol))
            @test isempty(_failing(i -> _fwf_padding_is_zero(fwf, i), 1:n))
        end
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
    # The rectangle `ext` as a GeoInterface polygon, for GeometryOps
    gi_rect(ext) = (((x1, x2), (y1, y2)) = (ext.X, ext.Y);
        GI.Polygon([GI.LinearRing([(x1, y1), (x1, y2), (x2, y2), (x2, y1), (x1, y1)])]))
    # Whether the interiors of floe i and the rectangle overlap, by GeometryOps, as the
    # docstring of interiors_intersect describes it
    function overlaps_go(fwf, i, ext)
        p, rect = Subzero._floe_poly(fwf, i), gi_rect(ext)
        return GO.intersects(p, rect) && !GO.touches(p, rect)
    end
    intersects(fwf, i, ext) = GOG.interiors_intersect(Subzero._floe_poly(fwf, i), ext)
    @testset "interiors_intersect matches GeometryOps with $FT" for FT in (Float64, Float32)
        floes, n_edge_cases, rects = _overlap_test_floes(FT)
        fwf = Subzero.FixedWidthFloes(floes)
        east = GI.extent(PeriodicBoundary(East, FT; x0 = 0, xf = 1e4, y0 = 0, yf = 1e4).poly)
        @test (east.X..., east.Y...) == (FT(1e4), FT(1.5e4), FT(-5e3), FT(1.5e4))
        @test [intersects(fwf, i, east) for i in 1:n_edge_cases] ==
            [false, false, false, true, true, true, true, false, false, true, true, false]
        # (rectangle, floe) pairs where the two differ. In Float32 they may differ by
        # round-off, but don't for these floes.
        mismatches = [(rect, i) for rect in [east; rects], i in eachindex(floes)
            if intersects(fwf, i, rect) != overlaps_go(fwf, i, rect)]
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
        fwf = Subzero.FixedWidthFloes(floes)
        expected(ext) = [overlaps_go(fwf, i, ext) ? 0x01 : 0x00 for i in eachindex(floes)]
        # rectangles for which the kernel's result differs
        @test isempty(_failing(ext -> marks(ext) == expected(ext), rects))
    end

    #=
    The geometry functions that GeometryOpsGPU takes from GeometryOps, compared with
    GeometryOps on random polygons and points chosen to hit the edge cases: points on edges
    and vertices, edges that touch, overlap, or are almost parallel. The kernel tests run
    them, through Subzero's methods for `ArrayPolygon`s (`area_poly`, ...), on every backend
    and compare with the CPU.
    =#
    #= Random polygons around the origin, so that most pairs overlap: star-shaped polygons
    with 3 to 9 vertices, and rectangles and L-shapes on a 250 m grid, which share edges and
    vertices with each other. Some are clockwise. =#
    function _geom_test_polys(::Type{FT}; n = 150, rng = Xoshiro(1)) where FT
        polys = Vector{Vector{Tuple{FT, FT}}}()
        while length(polys) < n
            if rand(rng) < 0.4
                nv = rand(rng, 3:9)
                θs = sort(2π .* rand(rng, nv))
                rs = rand(rng) < 0.5 ? fill(1e3, nv) : 1e3 .* (0.3 .+ rand(rng, nv))
                c = 500 .* (rand(rng, 2) .- 0.5)
                push!(polys, [(FT(c[1] + r * cos(θ)), FT(c[2] + r * sin(θ))) for (r, θ) in zip(rs, θs)])
            else
                g = FT.(250 .* sort(rand(rng, -4:4, 3)))
                h = FT.(250 .* sort(rand(rng, -4:4, 3)))
                (allunique(g) && allunique(h)) || continue
                push!(polys, rand(rng) < 0.5 ?
                    [(g[1], h[1]), (g[3], h[1]), (g[3], h[3]), (g[1], h[3])] :
                    [(g[1], h[1]), (g[3], h[1]), (g[3], h[2]), (g[2], h[2]), (g[2], h[3]), (g[1], h[3])])
            end
            rand(rng) < 0.3 && reverse!(polys[end])
        end
        return polys
    end
    # The polygons as closed apolys in a 3D array, the layout of FixedWidthFloes.poly
    function _array_polys(polys::Vector{Vector{Tuple{FT, FT}}}) where FT
        coords = zeros(FT, length(polys), maximum(length, polys) + 1, 2)
        apolys = map(enumerate(polys)) do (i, p)
            for (j, (x, y)) in enumerate([p; [p[1]]])
                coords[i, j, 1], coords[i, j, 2] = x, y
            end
            GOG.ArrayPolygon(coords, i, length(p) + 1)
        end
        return coords, apolys
    end
    _gi_poly(p) = GI.Polygon([GI.LinearRing([p; [p[1]]])])
    #= Test points for polygon p: its vertices, points along its edges, grid points, and
    random points. Each is also returned as Float64 and as a Float64 midpoint `0.5 * (a + b)`,
    as collisions compute them. =#
    function _geom_test_points(rng, p::Vector{Tuple{FT, FT}}, polys) where FT
        u = rand(rng)
        pt = if u < 0.3
            rand(rng, rand(rng, polys))
        elseif u < 0.5
            (FT(250 * rand(rng, -4:4)), FT(125 * rand(rng, -10:10)))
        elseif u < 0.7
            v, w = rand(rng, p), rand(rng, p)
            v .+ (rand(rng, 0:4) / 4) .* (w .- v)
        else
            (FT(3e3 * rand(rng) - 1.5e3), FT(3e3 * rand(rng) - 1.5e3))
        end
        return (pt, Float64.(pt), (0.5 * (pt[1] + pt[1]), 0.5 * (pt[2] + pt[2])))
    end

    #=
    Whether GOG's point tests give the same results as GeometryOps'. With Float32 points
    and polygons, GOG's test for a point on an edge (`isapprox` of two cross products) has
    the tolerance of Float32, where GeometryOps converts to Float64, so points close to an
    edge can be on it for GOG only. For those, `coveredby`, `within`, and the sign of the
    distance may differ. Returns whether the results are the same (or may differ), and if
    they may differ, the distance to the polygon's edges relative to the diagonal of its
    bounding box, else `nothing`. A point is on an edge for `isapprox` if it is closer to
    the edge than about `sqrt(eps(Float32))` times the edge's length.
    =#
    function _same_point_tests(q, r, gp, ::Type{FT}) where FT
        d, d_go = GOG.signed_distance(q, r, FT), GO.signed_distance(q, gp, FT)
        isequal(abs(d), abs(d_go)) || return false, nothing
        if q isa NTuple{2, Float32} && GOG._point_filled_curve_orientation(q, r) == GOG.point_on
            ext = GOG.extent(r)
            return true, abs(d) / hypot(ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1])
        end
        same = isequal(GOG.coveredby(q, r), GO.coveredby(q, gp)) &&
            isequal(GOG.within(q, r), GO.within(q, gp)) && isequal(d, d_go)
        return same, nothing
    end

    @testset "Area, centroid, and point tests with $FT" for FT in (Float64, Float32)
        polys = _geom_test_polys(FT)
        _, apolys = _array_polys(polys)
        gps = _gi_poly.(polys)
        @test all(i -> isequal(GOG.area(apolys[i], FT), GO.area(gps[i], FT)), eachindex(polys))
        @test all(i -> isequal(GOG.centroid(apolys[i], FT), GO.centroid(gps[i], FT)), eachindex(polys))
        rng = Xoshiro(2)
        orientations = GOG.PointOrientation[]
        same, near = true, Float64[]
        for _ in 1:5000
            i = rand(rng, eachindex(polys))
            r, gp = apolys[i], gps[i]
            for pt in _geom_test_points(rng, polys[i], polys)
                push!(orientations, GOG._point_filled_curve_orientation(pt, r))
                s, d = _same_point_tests(pt, r, gp, FT)
                same &= s
                isnothing(d) || push!(near, d)
            end
        end
        @test same
        @test all(<=(sqrt(eps(Float32))), near)
        # The points hit all three cases, including many points exactly on an edge
        @test count(==(GOG.point_on), orientations) > 1000
        @test count(==(GOG.point_in), orientations) > 1000
        @test count(==(GOG.point_out), orientations) > 1000
    end

    # All pairs of edges of the polygons, as (a1, a2, b1, b2): touching, crossing, and
    # collinear edges
    _polygon_edge_pairs(polys) = [(p[a], p[mod1(a + 1, end)], q[b], q[mod1(b + 1, end)])
        for p in polys for q in polys for a in eachindex(p) for b in eachindex(q)]
    #= Almost parallel edges that cross, as (a1, a2, b1, b2): rounding can put the computed
    crossing outside the edges' bounding boxes, and then the end closest to the other edge
    is used =#
    function _near_parallel_edge_pairs(::Type{FT}; n = 20000, rng = Xoshiro(5)) where FT
        return map(1:n) do _
            x0, y0 = 1e3 * rand(rng), 1e3 * rand(rng)
            θ, δ, L = 2π * rand(rng), 10.0^(-rand(rng, 4:12)), 1e3
            m = (x0 + L / 2 * cos(θ), y0 + L / 2 * sin(θ))
            a1, a2 = (FT(x0), FT(y0)), (FT(x0 + L * cos(θ)), FT(y0 + L * sin(θ)))
            b1 = (FT(m[1] - L / 3 * cos(θ + δ)), FT(m[2] - L / 3 * sin(θ + δ)))
            b2 = (FT(m[1] + L / 3 * cos(θ + δ)), FT(m[2] + L / 3 * sin(θ + δ)))
            (a1, a2, b1, b2)
        end
    end
    # The intersection of two edges, with its kind as an integer, from GeometryOps and GOG
    function _go_intersection_point(::Type{FT}, a1, a2, b1, b2) where FT
        kind, intr1, intr2 = GO._intersection_point(GO.Planar(), FT, (a1, a2), (b1, b2);
            exact = GO.True())
        return Int8(Int(kind)), intr1, intr2
    end
    function _gog_intersection_point(::Type{FT}, a1, a2, b1, b2) where FT
        kind, intr1, intr2 = GOG._intersection_point(FT, (a1, a2), (b1, b2))
        return Int8(kind), intr1, intr2
    end

    @testset "Edge intersections with $FT" for FT in (Float64, Float32)
        same_as_go(e) =
            isequal(_gog_intersection_point(FT, e...), _go_intersection_point(FT, e...))
        edges = _polygon_edge_pairs(_geom_test_polys(FT; n = 60))
        @test all(same_as_go, edges)
        kinds = [_gog_intersection_point(FT, e...)[1] for e in edges]
        @test all(k -> count(==(k), kinds) > 100,
            Int8.((GOG.line_cross, GOG.line_hinge, GOG.line_over, GOG.line_out)))
        edges = _near_parallel_edge_pairs(FT)
        @test all(same_as_go, edges)
        n_endpoints = count(edges) do e
            kind, (pt, _), _ = _gog_intersection_point(FT, e...)
            kind == Int8(GOG.line_cross) && pt in e
        end
        # In Float64, rounding hardly ever puts the crossing outside the bounding boxes
        FT == Float32 && @test n_endpoints > 100
    end

    #= With Float32 points and polygons, the maths is in Float32: no intermediate value in
    these functions (and the functions inlined into them) is a Float64 =#
    @testset "No Float64 maths with Float32" begin
        float64_values(f, types) = count(==(Float64), only(Base.code_typed(f, types))[1].ssavaluetypes)
        P = GOG.ArrayPolygon{Array{Float32, 3}}
        pt = NTuple{2, Float32}
        ext = typeof(Subzero.Extents.Extent(X = (1f0, 2f0), Y = (1f0, 2f0)))
        F = Type{Float32}
        for (f, types) in (
            (GOG.Predicates._orient, Tuple{GOG.True, pt, pt, pt}),
            (GOG.Predicates._orient, Tuple{GOG.False, pt, pt, pt}),
            (Subzero.GeometryOpsGPU.AdaptivePredicates.orient2, Tuple{pt, pt, pt}),
            (GOG._point_filled_curve_orientation, Tuple{pt, P}),
            (GOG.coveredby, Tuple{pt, P}),
            (GOG.within, Tuple{pt, P}),
            (GOG.area, Tuple{P, F}),
            (GOG.centroid, Tuple{P, F}),
            (GOG.signed_distance, Tuple{pt, P, F}),
            (GOG._intersection_point, Tuple{F, NTuple{2, pt}, NTuple{2, pt}}),
            (GOG._find_hinge_intersection, Tuple{F, pt, pt, pt, pt, Int, Int, Int}),
            (GOG._find_cross_intersection, Tuple{F, pt, pt, pt, pt, ext, ext}),
            (GOG._nearest_endpoint, Tuple{F, pt, pt, pt, pt}),
            (GOG.intersection_points!, Tuple{Array{Float32, 3}, Int, P, P, F}),
            (Subzero.get_intersection_points!, Tuple{Array{Float32, 3}, Int, P, P}),
        )
            @test float64_values(f, types) == 0
        end
    end

    @testset "Intersection points and matching vertices with $FT" for FT in (Float64, Float32)
        polys = _geom_test_polys(FT)
        _, apolys = _array_polys(polys)
        gps = _gi_poly.(polys)
        points = zeros(1, 64, 2)
        idxs = zeros(Int, 1, 64)
        same_points, same_idxs, n_nonempty = true, true, 0
        rng = Xoshiro(3)
        for _ in 1:3000
            i, j = rand(rng, eachindex(polys), 2)
            i == j && continue
            ref = GO.intersection_points(gps[i], gps[j])
            m = GOG.intersection_points!(points, 1, apolys[i], apolys[j], Float64)
            same_points &= isequal([(points[1, l, 1], points[1, l, 2]) for l in 1:m], ref)
            isempty(ref) && continue
            n_nonempty += 1
            c = Subzero.which_vertices_match_points!(idxs, 1, points, m, apolys[j])
            same_idxs &= isequal(idxs[1, 1:c], Subzero.which_vertices_match_points(ref, gps[j]))
        end
        @test same_points
        @test same_idxs
        @test n_nonempty > 1000
        #= Two squares that overlap in the region [1, 2] × [0, 1]: their edges cross at (1, 0)
        and (2, 1), which are vertices 1 and 3 of the region =#
        sq1 = [(FT(0), FT(0)), (FT(2), FT(0)), (FT(2), FT(2)), (FT(0), FT(2))]
        sq2 = [(FT(1), FT(-1)), (FT(3), FT(-1)), (FT(3), FT(1)), (FT(1), FT(1))]
        region = [(FT(1), FT(0)), (FT(2), FT(0)), (FT(2), FT(1)), (FT(1), FT(1))]
        _, (r1, r2, rr) = _array_polys([sq1, sq2, region])
        pts = zeros(1, 2, 2)
        @test GOG.intersection_points!(pts, 1, r1, r2, Float64) == 2
        @test pts[1, :, :] == [1 0; 2 1]
        idx = zeros(Int, 1, 2)
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 2, rr) == 2
        @test idx[1, :] == [1, 3]
        # Too little room for the points or indices: -1
        @test GOG.intersection_points!(zeros(1, 1, 2), 1, r1, r2, Float64) == -1
        @test Subzero.which_vertices_match_points!(zeros(Int, 1, 1), 1, pts, 2, rr) == -1
        # Like which_vertices_match_points, a single point matches nothing; no points either
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 1, rr) == 0
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 0, rr) == 0
    end

    #= Real floe shapes (7 to 591 vertices, from test/inputs/floe_shapes.jld2) and the
    synthetic polygons, at large coordinates, where floats are further apart: shifted by 1e6 m
    in x and y. Points include vertices and points on edges nudged by one or two floats, so
    the tolerance of the point-in-polygon test decides whether they are on an edge. Pairs
    overlap because one polygon is moved onto the other. =#
    @testset "Real floe shapes and large coordinates with $FT, offset $offset" for
            FT in (Float64, Float32), offset in (0.0, 1e6)
        real_shapes = jldopen(joinpath(@__DIR__, "inputs", "floe_shapes.jld2")) do f
            [[(c[1], c[2]) for c in p[1]] for p in vec(f["floe_vertices"])[1:4:end]]
        end
        base = [real_shapes; _geom_test_polys(Float64; n = 100)]
        shifted(p, Δx, Δy) = [(FT(x + Δx + offset), FT(y + Δy + offset)) for (x, y) in p]
        polys = [shifted(p, 0, 0) for p in base]
        _, apolys = _array_polys(polys)
        gps = _gi_poly.(polys)
        @test all(i -> isequal(GOG.area(apolys[i], FT), GO.area(gps[i], FT)), eachindex(polys))
        @test all(i -> isequal(GOG.centroid(apolys[i], FT), GO.centroid(gps[i], FT)), eachindex(polys))
        rng = Xoshiro(7)
        same, n_on, near = true, 0, Float64[]
        for (p, r, gp) in zip(polys, apolys, gps), _ in 1:10
            k = rand(rng, eachindex(p))
            e1, e2 = p[k], p[mod1(k + 1, end)]
            t = FT(rand(rng))
            on_edge = (e1[1] + t * (e2[1] - e1[1]), e1[2] + t * (e2[2] - e1[2]))
            v = rand(rng) < 0.5 ? e1 : on_edge
            pt = (nextfloat(v[1], rand(rng, -2:2)), nextfloat(v[2], rand(rng, -2:2)))
            for q in (pt, Float64.(pt))
                n_on += GOG._point_filled_curve_orientation(q, r) == GOG.point_on
                s, d = _same_point_tests(q, r, gp, FT)
                same &= s
                isnothing(d) || push!(near, d)
            end
        end
        @test same
        @test all(<=(sqrt(eps(Float32))), near)
        @test n_on > 100
        points, idxs = zeros(1, 64, 2), zeros(Int, 1, 64)
        same_points, same_idxs, n_nonempty = true, true, 0
        for i in eachindex(base), _ in 1:2
            j = rand(rng, eachindex(base))
            i == j && continue
            # move polygon j close to polygon i
            ci = (mean(first, base[i]), mean(last, base[i]))
            cj = (mean(first, base[j]), mean(last, base[j]))
            d = sqrt(abs(GO.area(_gi_poly(base[i])))) * (rand(rng) - 0.5)
            pj = shifted(base[j], ci[1] - cj[1] + d, ci[2] - cj[2] + rand(rng) * d)
            _, (ri, rj) = _array_polys([polys[i], pj])
            ref = GO.intersection_points(gps[i], _gi_poly(pj))
            m = GOG.intersection_points!(points, 1, ri, rj, Float64)
            same_points &= isequal([(points[1, l, 1], points[1, l, 2]) for l in 1:m], ref)
            isempty(ref) && continue
            n_nonempty += 1
            c = Subzero.which_vertices_match_points!(idxs, 1, points, m, rj)
            same_idxs &= isequal(idxs[1, 1:c], Subzero.which_vertices_match_points(ref, _gi_poly(pj)))
        end
        @test same_points
        @test same_idxs
        @test n_nonempty > 100
    end
    function _geom_point_tests!(orient, sdist, i, s, v, ::Type{FT}, point, r) where FT
        orient[i, s, v] = Int8(Subzero.GeometryOpsGPU._point_filled_curve_orientation(point, r)) +
            Int8(4) * Subzero.coveredby_poly(point, r) +
            Int8(8) * Subzero.GeometryOpsGPU.within(point, r)
        sdist[i, s, v] = Subzero.dist_to_poly(point, r, FT)
        return
    end
    #= Kernel test: each floe computes its area and centroid, the point tests for its test
    points, and the intersection points and matching vertices with a partner floe =#
    function _geom_kernel!(floes, i, pts, pts64, mids, orient, sdist, areas, cents, partner,
            ipts, nipts, idxs, nidxs)
        FT = eltype(floes.area)
        r = Subzero._floe_poly(floes, i)
        areas[i] = Subzero.area_poly(r, FT)
        cents[i, 1], cents[i, 2] = Subzero.centroid_poly(r, FT)
        # One call per point type: a loop over the three would be type-unstable for Float32
        for s in axes(pts, 2)
            _geom_point_tests!(orient, sdist, i, s, 1, FT, (pts[i, s, 1], pts[i, s, 2]), r)
            _geom_point_tests!(orient, sdist, i, s, 2, FT, (pts64[i, s, 1], pts64[i, s, 2]), r)
            _geom_point_tests!(orient, sdist, i, s, 3, FT, (mids[i, s, 1], mids[i, s, 2]), r)
        end
        rj = Subzero._floe_poly(floes, partner[i])
        m = Subzero.get_intersection_points!(ipts, i, r, rj, Float64)
        nipts[i] = m
        nidxs[i] = m < 0 ? 0 : Subzero.which_vertices_match_points!(idxs, i, ipts, m, rj)
        return
    end
    @testset "Kernels on $backend with $FT" for backend in test_backends(), FT in (Float64, Float32)
        polys = _geom_test_polys(FT; n = 100)
        fs = FloeSettings(FT; min_floe_area = 0)
        rng = Xoshiro(4)
        floes = StructArray([
            Floe{FT}([[[x, y] for (x, y) in [p; [p[1]]]]], FT(0.25); floe_settings = fs, rng)
            for p in polys
        ])
        fwf = Subzero.FixedWidthFloes(floes)
        n, ns = length(floes), 6
        pts, pts64, mids = zeros(FT, n, ns, 2), zeros(n, ns, 2), zeros(n, ns, 2)
        for i in 1:n, s in 1:ns
            p, p64, mid = _geom_test_points(rng, polys[i], polys)
            pts[i, s, :] .= p
            pts64[i, s, :] .= p64
            mids[i, s, :] .= mid
        end
        partner = Int32[mod1(i + rand(rng, 1:(n - 1)), n) for i in 1:n]
        outputs() = (zeros(Int8, n, ns, 3), zeros(FT, n, ns, 3), zeros(FT, n), zeros(FT, n, 2),
            zeros(n, 32, 2), zeros(Int32, n), zeros(Int32, n, 32), zeros(Int32, n))
        args(out) = (pts, pts64, mids, out[1], out[2], out[3], out[4], partner, out[5:8]...)
        # Reference: the same function called on the CPU without a kernel
        ref = outputs()
        for i in 1:n
            _geom_kernel!(fwf, i, args(ref)...)
        end
        dev = Subzero.adapt(backend, outputs())
        dev_args = Subzero.adapt(backend, args(dev))
        Subzero.launch_per_floe!(_geom_kernel!, backend, Subzero.adapt(backend, fwf), dev_args...)
        Subzero.KernelAbstractions.synchronize(backend)
        got = Subzero.adapt(Array, dev)
        @test all(isequal.(got, ref))
        # The partners do intersect, so the intersection points are tested
        @test count(>(0), ref[6]) > n ÷ 2
    end

    #= Kernel test of edge intersections. Its purpose is to make sure that GPUs give the same
    intersections as GeometryOps on every branch of `GOG._intersection_point`, including the rare
    ones. A GPU compiler could do floating-point maths differently from the CPU (fuse a
    multiplication and an addition, use a less exact square root, or reorder operations),
    which would change the results, and could break the exact orientation test, which
    relies on IEEE arithmetic. The per-floe kernel test above has too few edge pairs to
    reach the rare branches, such as the nearest end of almost parallel edges, so this test
    runs all the edge pairs of the CPU test above (about 150,000) in a kernel and compares
    them with GeometryOps directly. Kernels name the module in full, because the local `GOG`
    would be captured, and a module can't be passed to a kernel.

    Each thread computes the intersection of one pair: `launch_per_floe!` only needs
    `n_floes`, so it gets a NamedTuple with the number of pairs instead of floes. =#
    function _edge_intersection_kernel!(_, i, edges, kinds, intrs)
        FT = eltype(edges)
        a1, a2 = (edges[i, 1], edges[i, 2]), (edges[i, 3], edges[i, 4])
        b1, b2 = (edges[i, 5], edges[i, 6]), (edges[i, 7], edges[i, 8])
        kind, (p1, (α1, β1)), (p2, (α2, β2)) =
            Subzero.GeometryOpsGPU._intersection_point(FT, (a1, a2), (b1, b2))
        kinds[i] = Int8(kind)
        intrs[i, 1], intrs[i, 2], intrs[i, 3], intrs[i, 4] = p1[1], p1[2], α1, β1
        intrs[i, 5], intrs[i, 6], intrs[i, 7], intrs[i, 8] = p2[1], p2[2], α2, β2
        return
    end
    @testset "Edge intersections on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        edges = [_polygon_edge_pairs(_geom_test_polys(FT; n = 60));
            _near_parallel_edge_pairs(FT)]
        n = length(edges)
        edge_array = [e[k][d] for e in edges, k in 1:4, d in 1:2]
        edge_array = reshape(permutedims(edge_array, (1, 3, 2)), n, 8)
        dev_kinds = Subzero.KernelAbstractions.zeros(backend, Int8, n)
        dev_intrs = Subzero.KernelAbstractions.zeros(backend, FT, n, 8)
        Subzero.launch_per_floe!(_edge_intersection_kernel!, backend, (; n_floes = n),
            Subzero.adapt(backend, edge_array), dev_kinds, dev_intrs)
        Subzero.KernelAbstractions.synchronize(backend)
        kinds, intrs = Array(dev_kinds), Array(dev_intrs)
        got(i) = (kinds[i], ((intrs[i, 1], intrs[i, 2]), (intrs[i, 3], intrs[i, 4])),
            ((intrs[i, 5], intrs[i, 6]), (intrs[i, 7], intrs[i, 8])))
        @test all(i -> isequal(got(i), _go_intersection_point(FT, edges[i]...)), 1:n)
    end
end
