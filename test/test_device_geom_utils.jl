#=
The geometry functions in `src/tools/device_geom_utils.jl` replace GeometryOps functions in
GPU kernels, and must give the same results bit for bit. The first testsets compare them with
GeometryOps on the CPU, on random polygons and points chosen to hit the edge cases: points on
edges and vertices, edges that touch, overlap, or are almost parallel. The last testset runs
them in a kernel on every backend and compares with the CPU.
=#
@testset "Device geometry functions" begin
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
    # GeometryOps' intersection of two edges, with its kind as one of our _LINE_* values
    function _go_intersection_point(::Type{FT}, a1, a2, b1, b2) where FT
        kind, intr1, intr2 = GO._intersection_point(GO.Planar(), FT, (a1, a2), (b1, b2);
            exact = GO.True())
        return Int8(Int(kind)), intr1, intr2
    end

    @testset "Edge intersections with $FT" for FT in (Float64, Float32)
        same_as_go(e) =
            isequal(Subzero._intersection_point(FT, e...), _go_intersection_point(FT, e...))
        edges = _polygon_edge_pairs(_geom_test_polys(FT; n = 60))
        @test all(same_as_go, edges)
        kinds = [Subzero._intersection_point(FT, e...)[1] for e in edges]
        @test all(k -> count(==(k), kinds) > 100,
            (Subzero._LINE_CROSS, Subzero._LINE_HINGE, Subzero._LINE_OVER, Subzero._LINE_OUT))
        edges = _near_parallel_edge_pairs(FT)
        @test all(same_as_go, edges)
        n_endpoints = count(edges) do e
            kind, (pt, _), _ = Subzero._intersection_point(FT, e...)
            kind == Subzero._LINE_CROSS && pt in e
        end
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
        @test all(i -> isequal(Subzero.area_poly(apolys[i], FT), GO.area(gps[i], FT)), eachindex(polys))
        @test all(i -> isequal(Subzero.centroid_poly(apolys[i], FT), GO.centroid(gps[i], FT)), eachindex(polys))
        rng = Xoshiro(7)
        same, n_on = trues(3), 0
        for (p, r, gp) in zip(polys, apolys, gps), _ in 1:10
            k = rand(rng, eachindex(p))
            e1, e2 = p[k], p[mod1(k + 1, end)]
            t = FT(rand(rng))
            on_edge = (e1[1] + t * (e2[1] - e1[1]), e1[2] + t * (e2[2] - e1[2]))
            v = rand(rng) < 0.5 ? e1 : on_edge
            pt = (nextfloat(v[1], rand(rng, -2:2)), nextfloat(v[2], rand(rng, -2:2)))
            for q in (pt, Float64.(pt))
                n_on += Subzero._point_poly_orientation(q, r) == Subzero._POINT_ON
                same[1] &= isequal(Subzero.coveredby_poly(q, r), GO.coveredby(q, gp))
                same[2] &= isequal(Subzero._within_poly(q, r), GO.within(q, gp))
                same[3] &= isequal(Subzero.dist_to_poly(q, r, FT), GO.signed_distance(q, gp, FT))
            end
        end
        @test all(same)
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
            m = Subzero.get_intersection_points!(points, 1, ri, rj, Float64)
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

    #= Kernel test of edge intersections. Its purpose is to make sure that GPUs give the same
    intersections as GeometryOps on every branch of `_intersection_point`, including the rare
    ones. A GPU compiler could do floating-point maths differently from the CPU (fuse a
    multiplication and an addition, use a less exact square root, or reorder operations),
    which would change the results, and could break the exact orientation test, which
    relies on IEEE arithmetic. The per-floe kernel test above has too few edge pairs to
    reach the rare branches, such as the nearest end of almost parallel edges, so this test
    runs all the edge pairs of the CPU test above (about 150,000) in a kernel and compares
    them with GeometryOps directly.

    Each thread computes the intersection of one pair: `launch_per_floe!` only needs
    `n_floes`, so it gets a NamedTuple with the number of pairs instead of floes. =#
    function _edge_intersection_kernel!(_, i, edges, kinds, intrs)
        FT = eltype(edges)
        a1, a2 = (edges[i, 1], edges[i, 2]), (edges[i, 3], edges[i, 4])
        b1, b2 = (edges[i, 5], edges[i, 6]), (edges[i, 7], edges[i, 8])
        kind, (p1, (α1, β1)), (p2, (α2, β2)) =
            Subzero._intersection_point(FT, a1, a2, b1, b2)
        kinds[i] = kind
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

    #= Whether the device `check_intersects` gives the same as GeometryOps for polygons p
    and q (vectors of (x, y), not closed), and whether the pair has edges that only touch or
    overlap, which `check_intersects` handles by splitting edges =#
    function _compare_intersects(p, q)
        _, (pa, pb) = _array_polys([p, q])
        same = Subzero.check_intersects(pa, pb) == Subzero.check_intersects(_gi_poly(p), _gi_poly(q))
        T = Float64
        edge(r, i) = ((T(r[i][1]), T(r[i][2])), (T(r[mod1(i + 1, end)][1]), T(r[mod1(i + 1, end)][2])))
        touching = any(Iterators.product(eachindex(p), eachindex(q))) do (i, j)
            kind, _, _ = Subzero._intersection_point(T, edge(p, i)..., edge(q, j)...,
                Subzero._orient_inexact)
            kind == Subzero._LINE_HINGE || kind == Subzero._LINE_OVER
        end
        return same, touching
    end

    @testset "Polygons intersecting with $FT, offset $offset" for FT in (Float64, Float32),
            offset in (0.0, 1e6)
        square(x0, y0, s) = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s)]
        l_shape = [(0, 0), (3, 0), (3, 1), (1, 1), (1, 3), (0, 3)]
        special_pairs = (
            crossing = (square(0, 0, 2), square(1, 1, 2)),
            identical = (square(0, 0, 2), square(0, 0, 2)),
            shared_edge = (square(0, 0, 2), square(2, 0, 2)),
            partly_shared_edge = (square(0, 0, 2), square(1, -2, 2)),
            touching_corners = (square(0, 0, 2), square(2, 2, 2)),
            vertex_on_edge = (square(0, 0, 2), [(1, 2), (2, 3), (0, 3)]),
            nested = (square(0, 0, 4), square(1, 1, 2)),
            nested_sharing_edge = (square(0, 0, 4), square(0, 1, 2)),
            disjoint = (square(0, 0, 1), square(3, 3, 1)),
            # bounding boxes overlap, but the polygons don't touch
            in_notch = (l_shape, square(1.5, 1.5, 1)),
            touching_in_notch = (l_shape, square(1, 1, 1)),
        )
        rng = Xoshiro(12)
        test_pairs = [([Float64.(v) for v in p], [Float64.(v) for v in q]) for (p, q) in special_pairs]
        synthetic = _geom_test_polys(Float64; n = 60)
        append!(test_pairs, [(synthetic[i], synthetic[j]) for i in eachindex(synthetic)
            for j in i:length(synthetic)])
        # real floe shapes next to each other, close enough that some touch
        real_shapes = _floe_shape_rings()[1:8:end]
        for p in real_shapes
            q = rand(rng, real_shapes)
            cp = (mean(first, p), mean(last, p))
            cq = (mean(first, q), mean(last, q))
            Δx = cp[1] - cq[1] + (maximum(first, p) - minimum(first, q)) * (0.9 + 0.2 * rand(rng))
            push!(test_pairs, (p, [(x + Δx, y + cp[2] - cq[2]) for (x, y) in q]))
        end
        #= As in calc_normal_force: overlap regions of real floe shapes, and the regions after
        moving one floe a little, which share vertices and edges with the first ones. These
        are Float64, as intersect_polys returns. =#
        regions = Vector{Tuple{Float64, Float64}}[]
        for p in real_shapes[1:2:end]
            q = rand(rng, real_shapes)
            cp = (mean(first, p), mean(last, p))
            cq = (mean(first, q), mean(last, q))
            q = [(x + cp[1] - cq[1] + 300 * randn(rng), y + cp[2] - cq[2]) for (x, y) in q]
            region_polys(p, q) = [collect(GI.getpoint(GI.getexterior(r)))[1:(end - 1)]
                for r in Subzero.intersect_polys(_gi_poly(p), _gi_poly(q), Float64)]
            moved = [(x + 20 * randn(rng), y) for (x, y) in p]
            for r in region_polys(p, q), s in region_polys(moved, q)
                push!(regions, r, s)
            end
        end
        shift(p) = [(FT(x + offset), FT(y + offset)) for (x, y) in p]
        shift64(p) = [(x + offset, y + offset) for (x, y) in p]
        same, n_true, n_false, n_touching = true, 0, 0, 0
        all_pairs = [[(shift(p), shift(q)) for (p, q) in test_pairs];
            [(shift64(regions[i]), shift64(regions[i + 1])) for i in 1:2:length(regions)]]
        for (p, q) in all_pairs, (r, s) in ((p, q), (q, p))
            same_rs, touching = _compare_intersects(r, s)
            same &= same_rs
            result = Subzero.check_intersects(_gi_poly(r), _gi_poly(s))
            n_true += result
            n_false += !result
            n_touching += touching
        end
        @test same
        @test n_true > 1000
        @test n_false > 200
        @test n_touching > 300
    end

    @testset "Translating polygons with $FT" for FT in (Float64, Float32)
        polys = [_geom_test_polys(FT; n = 20); [[(FT(x), FT(y)) for (x, y) in p]
            for p in _floe_shape_rings()[1:16:end]]]
        coords, apolys = _array_polys(polys)
        rng = Xoshiro(13)
        dest = zeros(FT, 1, size(coords, 2), 2)
        same = true
        for (p, ap) in zip(polys, apolys), Δ in (FT.(100 .* randn(rng, 2)), 100 .* randn(rng, 2))
            moved = Subzero._ArrayPoly(dest, 1, Subzero._translate_poly!(dest, 1, ap, Δ...))
            expected = collect(GI.getpoint(GI.getexterior(
                Subzero._translate_poly(FT, _gi_poly(p), Δ...))))
            got = [Subzero._getpoint(moved, j) for j in 1:Subzero._npoint(moved)]
            same &= isequal(got, expected)
        end
        @test same
        @test Subzero._translate_poly!(zeros(FT, 1, 4, 2), 1, apolys[1], 1, 1) == -1
    end

    #= Kernel test: each floe checks whether it intersects a partner floe, and the partner
    moved by a small distance, as in calc_normal_force =#
    function _intersects_kernel!(floes, i, partner, shifts, moved, results)
        poly, partner_poly = Subzero._floe_poly(floes, i), Subzero._floe_poly(floes, partner[i])
        results[i, 1] = Subzero.check_intersects(poly, partner_poly)
        n = Subzero._translate_poly!(moved, i, partner_poly, shifts[i, 1], shifts[i, 2])
        moved_poly = Subzero._ArrayPoly(moved, i, n)
        results[i, 2] = Subzero.check_intersects(poly, moved_poly)
        return
    end
    @testset "Polygons intersecting on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        polys = _geom_test_polys(FT; n = 100)
        fs = FloeSettings(FT; min_floe_area = 0)
        rng = Xoshiro(14)
        floes = StructArray([
            Floe{FT}([[[x, y] for (x, y) in [p; [p[1]]]]], FT(0.25); floe_settings = fs, rng)
            for p in polys
        ])
        fwf = Subzero.FixedWidthFloes(floes)
        n = length(floes)
        partner = Int32[mod1(i + rand(rng, 1:(n - 1)), n) for i in 1:n]
        # some shifts move the partner far enough not to intersect, and the grid polygons by
        # whole grid cells, so that they touch
        shifts = FT.(rand(rng, (0, 250, 1000, 2000), n, 2) .* rand(rng, (-1, 1), n, 2))
        outputs() = (zeros(FT, n, size(fwf.poly, 2), 2), zeros(Bool, n, 2))
        ref = outputs()
        for i in 1:n
            _intersects_kernel!(fwf, i, partner, shifts, ref...)
        end
        moved(i) = Subzero._translate_poly(FT, floes.poly[partner[i]], shifts[i, :]...)
        expected = [Subzero.check_intersects(floes.poly[i],
            k == 1 ? floes.poly[partner[i]] : moved(i)) for i in 1:n, k in 1:2]
        @test ref[2] == expected
        dev = Subzero.adapt(backend, outputs())
        Subzero.launch_per_floe!(_intersects_kernel!, backend, Subzero.adapt(backend, fwf),
            Subzero.adapt(backend, partner), Subzero.adapt(backend, shifts), dev...)
        Subzero.KernelAbstractions.synchronize(backend)
        got = Subzero.adapt(Array, dev)
        @test got[2] == ref[2]
        @test isequal(got[1], ref[1])
        @test 0 < count(ref[2]) < 2n
    end
end
