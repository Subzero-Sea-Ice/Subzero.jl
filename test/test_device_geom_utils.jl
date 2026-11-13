#=
The geometry functions in `src/tools/device_geom_utils.jl` replace GeometryOps functions in
GPU kernels, and must give the same results bit for bit. The first testsets compare them with
GeometryOps on the CPU, on random polygons and points chosen to hit the edge cases: points on
edges and vertices, edges that touch, overlap, or are almost parallel. The last testset runs
them in a kernel on every backend and compares with the CPU.
=#
@testset "Device geometry functions" begin
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
            Subzero._ArrayPoly(coords, i, length(p) + 1)
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

    @testset "Area, centroid, and point tests with $FT" for FT in (Float64, Float32)
        polys = _geom_test_polys(FT)
        _, apolys = _array_polys(polys)
        gps = _gi_poly.(polys)
        @test all(i -> isequal(Subzero.area_poly(apolys[i], FT), GO.area(gps[i], FT)), eachindex(polys))
        @test all(i -> isequal(Subzero.centroid_poly(apolys[i], FT), GO.centroid(gps[i], FT)), eachindex(polys))
        rng = Xoshiro(2)
        orientations = Int8[]
        same = trues(4)
        for _ in 1:5000
            i = rand(rng, eachindex(polys))
            r, gp = apolys[i], gps[i]
            for pt in _geom_test_points(rng, polys[i], polys)
                push!(orientations, Subzero._point_poly_orientation(pt, r))
                same[1] &= isequal(Subzero.coveredby_poly(pt, r), GO.coveredby(pt, gp))
                same[2] &= isequal(Subzero._within_poly(pt, r), GO.within(pt, gp))
                same[3] &= isequal(Subzero.dist_to_poly(pt, r, FT), GO.signed_distance(pt, gp, FT))
                same[4] &= isequal(Subzero._distance_curve(FT, pt, r), abs(GO.signed_distance(pt, gp, FT)))
            end
        end
        @test all(same)
        # The points hit all three cases, including many points exactly on an edge
        @test count(==(Subzero._POINT_ON), orientations) > 1000
        @test count(==(Subzero._POINT_IN), orientations) > 1000
        @test count(==(Subzero._POINT_OUT), orientations) > 1000
    end

    @testset "Edge intersections with $FT" for FT in (Float64, Float32)
        go_intersection(a1, a2, b1, b2) =
            GO._intersection_point(GO.Planar(), FT, (a1, a2), (b1, b2); exact = GO.True())
        same_as_go(a1, a2, b1, b2) = (got = Subzero._intersection_point(FT, a1, a2, b1, b2);
            ref = go_intersection(a1, a2, b1, b2);
            isequal((Int(got[1]), got[2], got[3]), (Int(ref[1]), ref[2], ref[3])))
        # All pairs of edges of the test polygons: touching, crossing, and collinear edges
        polys = _geom_test_polys(FT; n = 60)
        kinds = Int8[]
        same = true
        for p in polys, q in polys, a in eachindex(p), b in eachindex(q)
            a1, a2, b1, b2 = p[a], p[mod1(a + 1, end)], q[b], q[mod1(b + 1, end)]
            push!(kinds, Subzero._intersection_point(FT, a1, a2, b1, b2)[1])
            same &= same_as_go(a1, a2, b1, b2)
        end
        @test same
        @test all(k -> count(==(k), kinds) > 100,
            (Subzero._LINE_CROSS, Subzero._LINE_HINGE, Subzero._LINE_OVER, Subzero._LINE_OUT))
        #= Almost parallel edges that cross: rounding can put the computed crossing outside
        the edges' bounding boxes, and then the end closest to the other edge is used =#
        rng = Xoshiro(5)
        same = true
        n_endpoints = 0
        for _ in 1:20000
            x0, y0 = 1e3 * rand(rng), 1e3 * rand(rng)
            θ, δ, L = 2π * rand(rng), 10.0^(-rand(rng, 4:12)), 1e3
            m = (x0 + L / 2 * cos(θ), y0 + L / 2 * sin(θ))
            a1, a2 = (FT(x0), FT(y0)), (FT(x0 + L * cos(θ)), FT(y0 + L * sin(θ)))
            b1 = (FT(m[1] - L / 3 * cos(θ + δ)), FT(m[2] - L / 3 * sin(θ + δ)))
            b2 = (FT(m[1] + L / 3 * cos(θ + δ)), FT(m[2] + L / 3 * sin(θ + δ)))
            kind, (pt, _), _ = Subzero._intersection_point(FT, a1, a2, b1, b2)
            n_endpoints += kind == Subzero._LINE_CROSS && pt in (a1, a2, b1, b2)
            same &= same_as_go(a1, a2, b1, b2)
        end
        @test same
        # In Float64, rounding hardly ever puts the crossing outside the bounding boxes
        FT == Float32 && @test n_endpoints > 100
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
            m = Subzero.get_intersection_points!(points, 1, apolys[i], apolys[j], Float64)
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
        @test Subzero.get_intersection_points!(pts, 1, r1, r2, Float64) == 2
        @test pts[1, :, :] == [1 0; 2 1]
        idx = zeros(Int, 1, 2)
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 2, rr) == 2
        @test idx[1, :] == [1, 3]
        # Too little room for the points or indices: -1
        @test Subzero.get_intersection_points!(zeros(1, 1, 2), 1, r1, r2, Float64) == -1
        @test Subzero.which_vertices_match_points!(zeros(Int, 1, 1), 1, pts, 2, rr) == -1
        # Like which_vertices_match_points, a single point matches nothing; no points either
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 1, rr) == 0
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 0, rr) == 0
    end

    function _geom_point_tests!(orient, sdist, i, s, v, ::Type{FT}, point, r) where FT
        orient[i, s, v] = Subzero._point_poly_orientation(point, r) +
            Int8(4) * Subzero.coveredby_poly(point, r) +
            Int8(8) * Subzero._within_poly(point, r)
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
end
