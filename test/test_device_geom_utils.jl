#=
Tests of Subzero's geometry functions for polygons stored in arrays
(`src/tools/device_geom_utils.jl`) and of `_floe_poly`. The functions of GeometryOpsGPU that
they call are tested in `test_geometry_ops_gpu*.jl`; here, the wrappers are compared with
GeometryOpsGPU, `which_vertices_match_points!` with its CPU version, and all of them are run
in a kernel on floes, as the ported physics will.
=#
@testset "Device geometry functions" begin
    @testset "_floe_poly with $FT" for FT in (Float64, Float32)
        floes = _make_timestep_test_floes(FT, FloeSettings(FT))
        fwf = Subzero.FixedWidthFloes(floes)
        n = length(floes)
        expected(i) = [Tuple(FT.(p)) for p in GI.getpoint(GI.getexterior(floes.poly[i]))]
        @test all(i -> _gog_points(Subzero._floe_poly(fwf, i)) == expected(i), 1:n)
        # it shares the floes' arrays: changing it changes the floe
        GOG.transform!(GOG.Translate(FT(1), FT(2)), Subzero._floe_poly(fwf, 1))
        @test fwf.poly[1, 1, 1] == expected(1)[1][1] + 1
        @test fwf.poly[1, 1, 2] == expected(1)[1][2] + 2
    end

    # The extent of a periodic boundary, which add_ghosts! tests the floes with
    @testset "Extent of a periodic boundary with $FT" for FT in (Float64, Float32)
        east = GI.extent(PeriodicBoundary(East, FT; x0 = 0, xf = 1e4, y0 = 0, yf = 1e4).poly)
        @test (east.X..., east.Y...) == (FT(1e4), FT(1.5e4), FT(-5e3), FT(1.5e4))
        _, (inside, crossing) = _array_polys([
            [(FT(1e3), FT(1e3)), (FT(3e3), FT(1e3)), (FT(3e3), FT(3e3))],
            [(FT(8e3), FT(1e3)), (FT(11e3), FT(1e3)), (FT(11e3), FT(3e3))],
        ])
        @test !GOG.interiors_intersect(inside, east)
        @test GOG.interiors_intersect(crossing, east)
    end

    #= which_vertices_match_points! for the intersection points of pairs of polygons, compared
    with which_vertices_match_points on the same points =#
    function _same_matching_vertices!(idxs, points, m, region, gi_region)
        ref = [(points[1, l, 1], points[1, l, 2]) for l in 1:m]
        c = Subzero.which_vertices_match_points!(idxs, 1, points, m, region)
        return isequal(idxs[1, 1:c], Subzero.which_vertices_match_points(ref, gi_region))
    end
    @testset "which_vertices_match_points! with $FT" for FT in (Float64, Float32)
        polys = _geom_test_polys(FT)
        _, apolys = _array_polys(polys)
        gps = _gi_poly.(polys)
        points, idxs = zeros(FT, 1, 64, 2), zeros(Int, 1, 64)
        same, n_nonempty = true, 0
        rng = Xoshiro(3)
        for _ in 1:3000
            i, j = rand(rng, eachindex(polys), 2)
            i == j && continue
            m = Subzero.get_intersection_points!(points, 1, apolys[i], apolys[j])
            m > 0 || continue
            n_nonempty += 1
            same &= _same_matching_vertices!(idxs, points, m, apolys[j], gps[j])
        end
        @test same
        @test n_nonempty > 1000
        #= Two squares that overlap in the region [1, 2] × [0, 1]: their edges cross at (1, 0)
        and (2, 1), which are vertices 1 and 3 of the region =#
        sq1 = [(FT(0), FT(0)), (FT(2), FT(0)), (FT(2), FT(2)), (FT(0), FT(2))]
        sq2 = [(FT(1), FT(-1)), (FT(3), FT(-1)), (FT(3), FT(1)), (FT(1), FT(1))]
        region = [(FT(1), FT(0)), (FT(2), FT(0)), (FT(2), FT(1)), (FT(1), FT(1))]
        _, (r1, r2, rr) = _array_polys([sq1, sq2, region])
        pts = zeros(FT, 1, 2, 2)
        @test Subzero.get_intersection_points!(pts, 1, r1, r2) == 2
        idx = zeros(Int, 1, 2)
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 2, rr) == 2
        @test idx[1, :] == [1, 3]
        # Too little room for the indices: -1
        @test Subzero.which_vertices_match_points!(zeros(Int, 1, 1), 1, pts, 2, rr) == -1
        # Like which_vertices_match_points, a single point matches nothing; no points either
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 1, rr) == 0
        @test Subzero.which_vertices_match_points!(idx, 1, pts, 0, rr) == 0
    end

    # The same with real floe shapes at large coordinates, moved onto each other
    @testset "which_vertices_match_points! with real floe shapes, $FT, offset $offset" for
            FT in (Float64, Float32), offset in (0.0, 1e6)
        base = [_floe_shape_rings()[1:4:end]; _geom_test_polys(Float64; n = 100)]
        shifted(p, Δx, Δy) = [(FT(x + Δx + offset), FT(y + Δy + offset)) for (x, y) in p]
        polys = [shifted(p, 0, 0) for p in base]
        rng = Xoshiro(7)
        points, idxs = zeros(FT, 1, 64, 2), zeros(Int, 1, 64)
        same, n_nonempty = true, 0
        for i in eachindex(base), _ in 1:2
            j = rand(rng, eachindex(base))
            i == j && continue
            ci = (mean(first, base[i]), mean(last, base[i]))
            cj = (mean(first, base[j]), mean(last, base[j]))
            d = sqrt(abs(GO.area(_gi_poly(base[i])))) * (rand(rng) - 0.5)
            pj = shifted(base[j], ci[1] - cj[1] + d, ci[2] - cj[2] + rand(rng) * d)
            _, (ri, rj) = _array_polys([polys[i], pj])
            m = Subzero.get_intersection_points!(points, 1, ri, rj)
            m > 0 || continue
            n_nonempty += 1
            same &= _same_matching_vertices!(idxs, points, m, rj, _gi_poly(pj))
        end
        @test same
        @test n_nonempty > 100
    end

    # With Float32 floes, the maths is in Float32 (see "No Float64 maths with Float32" in
    # test_geometry_ops_gpu.jl)
    @testset "No Float64 maths with Float32" begin
        P = GOG.ArrayPolygon{Array{Float32, 3}}
        pt = NTuple{2, Float32}
        L = typeof(GOG.ClipLists(GOG.ClipLists(Float32, 1, 8), 1))
        Rg = typeof(GOG.PolyRegions(GOG.PolyRegions(Float32, 1, 1, 8), 1))
        F = Type{Float32}
        for (f, types) in (
            (Subzero.area_poly, Tuple{P, F}),
            (Subzero.centroid_poly, Tuple{P, F}),
            (Subzero.dist_to_poly, Tuple{pt, P, F}),
            (Subzero.coveredby_poly, Tuple{pt, P}),
            (Subzero.get_intersection_points!, Tuple{Array{Float32, 3}, Int, P, P}),
            (Subzero.intersect_polys!, Tuple{Rg, L, P, P}),
            (Subzero.which_vertices_match_points!, Tuple{Matrix{Int}, Int, Array{Float32, 3}, Int, P}),
        )
            @test _float64_values(f, types) == 0
        end
    end

    #= Kernel test: each floe calls the wrappers on itself, at its test points, and with a
    partner floe, as the ported physics will. intersect_polys! has a kernel test of its own
    below. =#
    function _wrappers_kernel!(floes, i, pts, out_f, out_c, ipts, idxs, partner)
        FT = eltype(floes.area)
        r, rj = Subzero._floe_poly(floes, i), Subzero._floe_poly(floes, partner[i])
        out_f[i, 1] = Subzero.area_poly(r, FT)
        out_f[i, 2], out_f[i, 3] = Subzero.centroid_poly(r, FT)
        for s in axes(pts, 2)
            pt = (pts[i, s, 1], pts[i, s, 2])
            out_f[i, 3 + s] = Subzero.dist_to_poly(pt, r, FT)
            out_c[i, s] = Subzero.coveredby_poly(pt, r)
        end
        m = Subzero.get_intersection_points!(ipts, i, r, rj)
        k = size(pts, 2)
        out_c[i, k + 1] = m
        out_c[i, k + 2] = m < 0 ? 0 : Subzero.which_vertices_match_points!(idxs, i, ipts, m, rj)
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
        n, k = length(floes), 6
        pts = zeros(FT, n, k, 2)
        for i in 1:n, s in 1:k
            pts[i, s, :] .= _geom_test_points(rng, polys[i], polys)
        end
        partner = [mod1(i + rand(rng, 1:(n - 1)), n) for i in 1:n]
        outputs() = (zeros(FT, n, 3 + k), zeros(Int, n, k + 2), zeros(FT, n, 32, 2),
            zeros(Int, n, 32))
        # Reference: the same function called on the CPU without a kernel
        ref = outputs()
        for i in 1:n
            _wrappers_kernel!(fwf, i, pts, ref..., partner)
        end
        # The wrappers give what GeometryOpsGPU gives
        r(i) = Subzero._floe_poly(fwf, i)
        @test all(i -> ref[1][i, 1] == GOG.area(r(i), FT), 1:n)
        @test all(i -> (ref[1][i, 2], ref[1][i, 3]) == GOG.centroid(r(i), FT), 1:n)
        @test all(i -> ref[1][i, 4] == GOG.signed_distance((pts[i, 1, 1], pts[i, 1, 2]), r(i), FT), 1:n)
        @test all(>=(0), ref[2][:, (k + 1):end])  # the buffers are large enough
        # intersect_polys! gives what GOG.intersection! gives
        clip(f, i) = (out = GOG.PolyRegions(GOG.PolyRegions(FT, 1, 8, 64), 1);
            m = f(out, GOG.ClipLists(GOG.ClipLists(FT, 1, 64), 1), r(i), r(partner[i]));
            (m, [_gog_points(GOG.region(out, l)) for l in 1:max(m, 0)]))
        @test all(i -> isequal(clip(Subzero.intersect_polys!, i), clip(GOG.intersection!, i)), 1:n)
        @test count(i -> clip(GOG.intersection!, i)[1] > 0, 1:n) > n ÷ 2
        dev = Subzero.adapt(backend, outputs())
        Subzero.launch_per_floe!(_wrappers_kernel!, backend, Subzero.adapt(backend, fwf),
            Subzero.adapt(backend, pts), dev..., Subzero.adapt(backend, partner))
        Subzero.KernelAbstractions.synchronize(backend)
        got = Subzero.adapt(Array, dev)
        @test all(isequal.(got, ref))
        # The partners do intersect, so the intersections are tested
        @test count(>(0), ref[2][:, k + 1]) > n ÷ 2
    end

    #= Kernel test of intersect_polys!: each floe intersects itself with a partner floe. Only
    with Float64 on JLBackend, to check that the wrapper works in a kernel: compiling the
    clipping into another kernel takes about 45 s, and the clipping itself
    (GOG.intersection!) is run in kernels on every backend, with both float types, in
    test_geometry_ops_gpu_clipping.jl. =#
    function _intersect_polys_kernel!(floes, i, lists, regions, counts, partner)
        counts[i] = Subzero.intersect_polys!(GOG.PolyRegions(regions, i), GOG.ClipLists(lists, i),
            Subzero._floe_poly(floes, i), Subzero._floe_poly(floes, partner[i]))
        return
    end
    @testset "intersect_polys! on JLBackend with Float64" begin
        backend, FT = JLArrays.JLBackend(), Float64
        polys = _geom_test_polys(FT; n = 100)
        fs = FloeSettings(FT; min_floe_area = 0)
        rng = Xoshiro(4)
        floes = StructArray([
            Floe{FT}([[[x, y] for (x, y) in [p; [p[1]]]]], FT(0.25); floe_settings = fs, rng)
            for p in polys
        ])
        fwf = Subzero.FixedWidthFloes(floes)
        n = length(floes)
        partner = [mod1(i + rand(rng, 1:(n - 1)), n) for i in 1:n]
        outputs() = (GOG.ClipLists(FT, n, 64), GOG.PolyRegions(FT, n, 8, 64), zeros(Int, n))
        # Reference: the same function called on the CPU without a kernel
        ref = outputs()
        for i in 1:n
            _intersect_polys_kernel!(fwf, i, ref..., partner)
        end
        dev = Subzero.adapt(backend, outputs())
        Subzero.launch_per_floe!(_intersect_polys_kernel!, backend, Subzero.adapt(backend, fwf),
            dev..., Subzero.adapt(backend, partner))
        Subzero.KernelAbstractions.synchronize(backend)
        got = Subzero.adapt(Array, dev)
        @test all(>=(0), ref[3])  # the buffers are large enough
        @test isequal(got[3], ref[3])
        regions(out, i, m) = [_gog_points(GOG.region(GOG.PolyRegions(out, i), l)) for l in 1:m]
        @test all(i -> isequal(regions(got[2], i, got[3][i]), regions(ref[2], i, ref[3][i])), 1:n)
        # most pairs overlap
        @test count(>(0), ref[3]) > n ÷ 2
    end

end
