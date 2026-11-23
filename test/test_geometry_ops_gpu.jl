#=
Tests of the GeometryOpsGPU submodule against GeometryOps, which accepts `ArrayPolygon`s as
GeoInterface polygons. With Float64 polygons the results must be the same, bit for bit; with
Float32 polygons they may differ by round-off. The tests use only GeometryOpsGPU, GeometryOps
and the packages they build on, with the helpers in `geometry_ops_gpu_utils.jl`, so that they
can move with GeometryOpsGPU.
=#
@testset "GeometryOpsGPU" begin
    Extent = Extents.Extent
    CT = CoordinateTransformations
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
        @test GI.extent(p) == GOG.extent(p) == Extent(X = (1.0, 2.0), Y = (2.0, 3.0))
    end

    # transform! changes the polygon in place and leaves the padding after its points alone
    @testset "transform! with $FT" for FT in (Float64, Float32)
        polys = _geom_test_polys(FT; n = 40)
        gps = _gi_poly.(polys)
        go_points(g) = [(FT(GI.x(q)), FT(GI.y(q))) for q in GI.getpoint(GI.getexterior(g))]
        # Translate adds Δx and Δy to every coordinate, as GO.transform with a Translation
        Δx, Δy = FT(-2e4), FT(1234.5)
        _, apolys = _array_polys(polys)
        foreach(p -> GOG.transform!(GOG.Translate(Δx, Δy), p), apolys)
        @test all(eachindex(polys)) do i
            _gog_points(apolys[i]) == go_points(GO.transform(CT.Translation(Δx, Δy), gps[i]))
        end
        @test all(_gog_padding_is_zero, apolys)
        #= RotateTranslate, as GO.transform with Translation ∘ recenter(LinearMap(Angle2d)): bit
        for bit in Float64, to round-off in Float32. The centres are tuples, as in _move_poly:
        with a Vector, StaticArrays multiplies with muladd, which may or may not be fused,
        depending on the CPU. =#
        Δα = FT(0.7)
        _, apolys = _array_polys(polys)
        centres = [GOG.centroid(p, FT) for p in apolys]
        foreach(i -> GOG.transform!(GOG.RotateTranslate(Δα, centres[i]..., Δx, Δy), apolys[i]),
            eachindex(apolys))
        expected(i) = go_points(GO.transform(CT.Translation(Δx, Δy) ∘
            CT.recenter(CT.LinearMap(Rotations.Angle2d(Δα)), centres[i]), gps[i]))
        @test all(eachindex(polys)) do i
            isapprox(reinterpret(FT, _gog_points(apolys[i])), reinterpret(FT, expected(i));
                rtol = FT == Float64 ? 0.0 : 10eps(FT))
        end
        @test all(_gog_padding_is_zero, apolys)
    end

    #=
    Polygons for testing `interiors_intersect`, and rectangles to test them with. Besides
    random star-shaped polygons, there are axis-aligned rectangles and L-shapes on a grid with
    a spacing of 250 m, and the rectangles' sides are on the same grid or go through polygon
    vertices, so polygons and rectangles often share edges or vertices. The polygons from
    `edge_cases` are tested with the rectangle 10 ≤ x ≤ 15 km, -5 ≤ y ≤ 15 km (in Subzero, the
    eastern periodic boundary of a 10 km × 10 km domain). Returns the polygons as vectors of
    points and as `ArrayPolygon`s, the number of edge cases, and the rectangles.
    =#
    function _overlap_test_polys(::Type{FT}; n = 200) where FT
        rng = Xoshiro(3)
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
        vertices = [p for poly in polys[(length(edge_cases) + 1):end] for p in poly]
        pick(i) = (u = rand(rng); u < 0.4 ? rand(rng, vertices)[i] :
            u < 0.7 ? FT(250 * rand(rng, -4:4)) : FT(4e3 * rand(rng) - 2e3))
        rects = Extent[]
        while length(rects) < 20
            xmin, xmax = minmax(pick(1), pick(1))
            ymin, ymax = minmax(pick(2), pick(2))
            xmin < xmax && ymin < ymax && push!(rects, Extent(X = (xmin, xmax), Y = (ymin, ymax)))
        end
        _, apolys = _array_polys(polys)
        return polys, apolys, length(edge_cases), rects
    end
    #= Whether the interiors of the polygon (a vector of points) and the rectangle overlap, by
    GeometryOps, as the docstring of interiors_intersect describes it =#
    function overlaps_go(p, ext, FT)
        gp, rect = _gi_poly(p), _gi_rect(ext)
        return GO.intersects(gp, rect) && !GO.touches(gp, rect)
    end
    @testset "interiors_intersect matches GeometryOps with $FT" for FT in (Float64, Float32)
        polys, apolys, n_edge_cases, rects = _overlap_test_polys(FT)
        east = Extent(X = (FT(1e4), FT(1.5e4)), Y = (FT(-5e3), FT(1.5e4)))
        @test [GOG.interiors_intersect(apolys[i], east) for i in 1:n_edge_cases] ==
            [false, false, false, true, true, true, true, false, false, true, true, false]
        # (rectangle, polygon) pairs where the two differ. In Float32 they may differ by
        # round-off, but don't for these polygons.
        mismatches = [(rect, i) for rect in [east; rects], i in eachindex(polys)
            if GOG.interiors_intersect(apolys[i], rect) != overlaps_go(polys[i], rect, FT)]
        @test isempty(mismatches)
    end
    # Bounds of another type than the polygon's aren't accepted
    @testset "interiors_intersect with Float32 polygons" begin
        _, apolys, _, _ = _overlap_test_polys(Float32)
        @test_throws MethodError GOG.interiors_intersect(apolys[1],
            Extent(X = (1e4, 1.5e4), Y = (-5e3, 1.5e4)))
    end
    # Kernel test: whether each polygon overlaps the rectangle
    _mark_overlap!(i, coords, ns, ext, marks) =
        (marks[i] = GOG.interiors_intersect(GOG.ArrayPolygon(coords, i, ns[i]), ext); nothing)
    @testset "interiors_intersect on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        polys, apolys, _, rects = _overlap_test_polys(FT)
        n = length(polys)
        coords = Adapt.adapt(backend, apolys[1].coords)
        ns = Adapt.adapt(backend, GOG.npoint.(apolys))
        function marks(ext)
            dev = KA.zeros(backend, Bool, n)
            _gog_launch!(_mark_overlap!, backend, n, coords, ns, ext, dev)
            return Array(dev)
        end
        expected(ext) = [overlaps_go(p, ext, FT) for p in polys]
        # rectangles for which the kernel's result differs
        @test isempty(filter(ext -> marks(ext) != expected(ext), rects))
    end

    #=
    The geometry functions that GeometryOpsGPU takes from GeometryOps, compared with
    GeometryOps on random polygons and points chosen to hit the edge cases: points on edges
    and vertices, edges that touch, overlap, or are almost parallel. The kernel tests run
    them on every backend and compare with the CPU.
    =#
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
        if FT == Float32 && GOG._point_filled_curve_orientation(q, r) == GOG.point_on
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
        for _ in 1:15000
            i = rand(rng, eachindex(polys))
            r, gp = apolys[i], gps[i]
            pt = _geom_test_points(rng, polys[i], polys)
            push!(orientations, GOG._point_filled_curve_orientation(pt, r))
            s, d = _same_point_tests(pt, r, gp, FT)
            same &= s
            isnothing(d) || push!(near, d)
        end
        @test same
        @test all(<=(sqrt(eps(Float32))), near)
        # Points of another type than the polygon's aren't accepted
        other = FT == Float64 ? Float32 : Float64
        pt = other.(polys[1][1])
        @test_throws MethodError GOG.coveredby(pt, apolys[1])
        @test_throws MethodError GOG.signed_distance(pt, apolys[1], FT)
        @test_throws MethodError GOG.area(apolys[1], other)
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
    #= The intersection of two edges, with its kind as an integer, from GeometryOps, with exact
    (default) or inexact orientation tests =#
    function _go_intersection_point(::Type{FT}, a1, a2, b1, b2; exact = true) where FT
        kind, intr1, intr2 = GO._intersection_point(GO.Planar(), FT, (a1, a2), (b1, b2);
            exact = exact ? GO.True() : GO.False())
        return Int8(Int(kind)), intr1, intr2
    end

    @testset "Intersection points with $FT" for FT in (Float64, Float32)
        polys = _geom_test_polys(FT)
        _, apolys = _array_polys(polys)
        gps = _gi_poly.(polys)
        points = zeros(FT, 1, 64, 2)
        same_points, n_nonempty = true, 0
        rng = Xoshiro(3)
        for _ in 1:3000
            i, j = rand(rng, eachindex(polys), 2)
            i == j && continue
            ref = GO.intersection_points(gps[i], gps[j], FT)
            m = GOG.intersection_points!(points, 1, apolys[i], apolys[j], FT)
            same_points &= isequal([(points[1, l, 1], points[1, l, 2]) for l in 1:m], ref)
            n_nonempty += !isempty(ref)
        end
        @test same_points
        @test n_nonempty > 1000
        # Two squares whose edges cross at (1, 0) and (2, 1)
        sq1 = [(FT(0), FT(0)), (FT(2), FT(0)), (FT(2), FT(2)), (FT(0), FT(2))]
        sq2 = [(FT(1), FT(-1)), (FT(3), FT(-1)), (FT(3), FT(1)), (FT(1), FT(1))]
        _, (r1, r2) = _array_polys([sq1, sq2])
        pts = zeros(FT, 1, 2, 2)
        @test GOG.intersection_points!(pts, 1, r1, r2, FT) == 2
        @test pts[1, :, :] == [1 0; 2 1]
        # Too little room for the points: -1
        @test GOG.intersection_points!(zeros(FT, 1, 1, 2), 1, r1, r2, FT) == -1
    end

    #= Real floe shapes (7 to 591 vertices, from test/inputs/floe_shapes.jld2) and the
    synthetic polygons, at large coordinates, where floats are further apart: shifted by 1e6 m
    in x and y. Points include vertices and points on edges nudged by one or two floats, so
    the tolerance of the point-in-polygon test decides whether they are on an edge. Pairs
    overlap because one polygon is moved onto the other. =#
    @testset "Real floe shapes and large coordinates with $FT, offset $offset" for
            FT in (Float64, Float32), offset in (0.0, 1e6)
        real_shapes = _floe_shape_rings()[1:4:end]
        base = [real_shapes; _geom_test_polys(Float64; n = 100)]
        shifted(p, Δx, Δy) = [(FT(x + Δx + offset), FT(y + Δy + offset)) for (x, y) in p]
        polys = [shifted(p, 0, 0) for p in base]
        _, apolys = _array_polys(polys)
        gps = _gi_poly.(polys)
        @test all(i -> isequal(GOG.area(apolys[i], FT), GO.area(gps[i], FT)), eachindex(polys))
        @test all(i -> isequal(GOG.centroid(apolys[i], FT), GO.centroid(gps[i], FT)), eachindex(polys))
        rng = Xoshiro(7)
        same, n_on, near = true, 0, Float64[]
        for (p, r, gp) in zip(polys, apolys, gps), _ in 1:20
            k = rand(rng, eachindex(p))
            e1, e2 = p[k], p[mod1(k + 1, end)]
            t = FT(rand(rng))
            on_edge = (e1[1] + t * (e2[1] - e1[1]), e1[2] + t * (e2[2] - e1[2]))
            v = rand(rng) < 0.5 ? e1 : on_edge
            pt = (nextfloat(v[1], rand(rng, -2:2)), nextfloat(v[2], rand(rng, -2:2)))
            n_on += GOG._point_filled_curve_orientation(pt, r) == GOG.point_on
            s, d = _same_point_tests(pt, r, gp, FT)
            same &= s
            isnothing(d) || push!(near, d)
        end
        @test same
        @test all(<=(sqrt(eps(Float32))), near)
        @test n_on > 100
        points = zeros(FT, 1, 64, 2)
        same_points, n_nonempty = true, 0
        for i in eachindex(base), _ in 1:2
            j = rand(rng, eachindex(base))
            i == j && continue
            # move polygon j close to polygon i
            ci = (mean(first, base[i]), mean(last, base[i]))
            cj = (mean(first, base[j]), mean(last, base[j]))
            d = sqrt(abs(GO.area(_gi_poly(base[i])))) * (rand(rng) - 0.5)
            pj = shifted(base[j], ci[1] - cj[1] + d, ci[2] - cj[2] + rand(rng) * d)
            _, (ri, rj) = _array_polys([polys[i], pj])
            ref = GO.intersection_points(gps[i], _gi_poly(pj), FT)
            m = GOG.intersection_points!(points, 1, ri, rj, FT)
            same_points &= isequal([(points[1, l, 1], points[1, l, 2]) for l in 1:m], ref)
            n_nonempty += !isempty(ref)
        end
        @test same_points
        @test n_nonempty > 100
    end

    #= The kernels below are launched with `_launch_float32_checked`, which checks that with
    Float32 they don't compute in Float64. Together they call every function of the module
    that kernels use, except the clipping (see `test_geometry_ops_gpu_clipping.jl`) and
    `interiors_intersect` and `transform!` (see `test_kernel_code.jl`). =#

    #= Kernel test: each polygon computes its area and centroid, the point tests for its test
    points, and the intersection points with a partner polygon =#
    function _geom_kernel!(i, coords, ns, pts, orient, sdist, areas, cents, partner, ipts, nipts)
        FT = eltype(coords)
        r = GOG.ArrayPolygon(coords, i, ns[i])
        areas[i] = GOG.area(r, FT)
        cents[i, 1], cents[i, 2] = GOG.centroid(r, FT)
        for s in axes(pts, 2)
            pt = (pts[i, s, 1], pts[i, s, 2])
            orient[i, s] = Int8(GOG._point_filled_curve_orientation(pt, r)) +
                Int8(4) * GOG.coveredby(pt, r) + Int8(8) * GOG.within(pt, r)
            sdist[i, s] = GOG.signed_distance(pt, r, FT)
        end
        j = partner[i]
        nipts[i] = GOG.intersection_points!(ipts, i, r, GOG.ArrayPolygon(coords, j, ns[j]), FT)
        return
    end
    @testset "Kernels on $backend with $FT" for backend in test_backends(), FT in (Float64, Float32)
        polys = _geom_test_polys(FT; n = 100)
        coords, apolys = _array_polys(polys)
        ns = GOG.npoint.(apolys)
        rng = Xoshiro(4)
        n, k = length(polys), 18
        pts = zeros(FT, n, k, 2)
        for i in 1:n, s in 1:k
            pts[i, s, :] .= _geom_test_points(rng, polys[i], polys)
        end
        partner = [mod1(i + rand(rng, 1:(n - 1)), n) for i in 1:n]
        outputs() = (zeros(Int8, n, k), zeros(FT, n, k), zeros(FT, n), zeros(FT, n, 2),
            zeros(FT, n, 32, 2), zeros(Int, n))
        args(out) = (coords, ns, pts, out[1], out[2], out[3], out[4], partner, out[5], out[6])
        # Reference: the same function called on the CPU without a kernel
        ref = outputs()
        for i in 1:n
            _geom_kernel!(i, args(ref)...)
        end
        dev = Adapt.adapt(backend, outputs())
        _launch_float32_checked(backend, FT) do
            _gog_launch!(_geom_kernel!, backend, n, Adapt.adapt(backend, args(dev))...)
        end
        got = Adapt.adapt(Array, dev)
        @test all(isequal.(got, ref))
        # The partners do intersect, so the intersection points are tested
        @test count(>(0), ref[6]) > n ÷ 2
    end

    #= Kernel test of edge intersections. Its purpose is to make sure that GPUs give the same
    intersections as GeometryOps on every branch of `GOG._intersection_point`, including the
    rare ones. A GPU compiler could do floating-point maths differently from the CPU (fuse a
    multiplication and an addition, use a less exact square root, or reorder operations),
    which would change the results, and could break the exact orientation test, which
    relies on IEEE arithmetic. The per-floe kernel test above has too few edge pairs to
    reach the rare branches, such as the nearest end of almost parallel edges, so this test
    runs all pairs of edges of 60 test polygons, which touch, cross, and overlap, and almost
    parallel edges (about 150,000 pairs) in a kernel and compares them with GeometryOps
    directly. Each thread computes the intersection of one pair, with exact or (as disjoint
    uses them) inexact orientation tests. =#
    function _edge_intersection_kernel!(i, edges, kinds, intrs, exact)
        FT = eltype(edges)
        a1, a2 = (edges[i, 1], edges[i, 2]), (edges[i, 3], edges[i, 4])
        b1, b2 = (edges[i, 5], edges[i, 6]), (edges[i, 7], edges[i, 8])
        kind, (p1, (α1, β1)), (p2, (α2, β2)) =
            GOG._intersection_point(FT, (a1, a2), (b1, b2); exact)
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
        # The intersections from the kernel, with exact or inexact orientation tests
        function kernel_results(exact)
            dev_kinds = KA.zeros(backend, Int8, n)
            dev_intrs = KA.zeros(backend, FT, n, 8)
            _launch_float32_checked(backend, FT) do
                _gog_launch!(_edge_intersection_kernel!, backend, n,
                    Adapt.adapt(backend, edge_array), dev_kinds, dev_intrs,
                    exact ? GOG.True() : GOG.False())
            end
            kinds, intrs = Array(dev_kinds), Array(dev_intrs)
            return [(kinds[i], ((intrs[i, 1], intrs[i, 2]), (intrs[i, 3], intrs[i, 4])),
                ((intrs[i, 5], intrs[i, 6]), (intrs[i, 7], intrs[i, 8]))) for i in 1:n]
        end
        ref = [_go_intersection_point(FT, e...) for e in edges]
        @test isequal(kernel_results(true), ref)
        #= With inexact orientation tests and Float32 edges, a point counts as on the other
        edge's line within Float32's tolerance (`isapprox`), where GeometryOps computes in
        Float64 (see `_same_point_tests`), so almost parallel edges can be classed
        differently, but only if GOG puts a point on the other edge's line =#
        inexact = kernel_results(false)
        on_line(e) = ((a1, a2, b1, b2) = e; any(iszero, (
            GOG.Predicates.orient(b1, b2, a1; exact = GOG.False()),
            GOG.Predicates.orient(b1, b2, a2; exact = GOG.False()),
            GOG.Predicates.orient(a1, a2, b1; exact = GOG.False()),
            GOG.Predicates.orient(a1, a2, b2; exact = GOG.False()))))
        differ = filter(i -> !isequal(inexact[i],
            _go_intersection_point(FT, edges[i]...; exact = false)), 1:n)
        @test FT == Float64 ? isempty(differ) : all(i -> on_line(edges[i]), differ)
        # The edge pairs reach every kind of intersection, and, in Float32, crossings that
        # rounding puts outside the edges' bounding boxes, for which the nearest end is used
        # (in Float64, that hardly ever happens)
        @test all(k -> count(r -> r[1] == Int8(k), ref) > 100,
            (GOG.line_cross, GOG.line_hinge, GOG.line_over, GOG.line_out))
        n_endpoints = count(i -> ref[i][1] == Int8(GOG.line_cross) && ref[i][2][1] in edges[i], 1:n)
        FT == Float32 && @test n_endpoints > 100
    end
end
