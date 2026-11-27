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

    @testset "_move_floe! on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        floes = _make_timestep_test_floes(FT, FloeSettings(FT))
        n = length(floes)
        Δx, Δy, Δα = FT(-31.4), FT(27.2), FT(0.7)
        fwf = _launch_on(Subzero._move_floe!, backend, Subzero.FixedWidthFloes(floes), Δx, Δy, Δα)
        expected = [Subzero._move_poly(FT, floes.poly[i], Δx, Δy, Δα, floes.centroid[i]...)
            for i in 1:n]
        @test isempty(_not_approx(
            [reinterpret(FT, _fwf_points(fwf, i)) for i in 1:n],
            [reinterpret(FT, collect(GI.getpoint(p))) for p in expected];
            rtol = FT == Float64 ? _float64_rtol(backend) : 10eps(FT),
        ))
        @test isempty(_not_approx([fwf.centroid[i, :] for i in 1:n],
            [floes.centroid[i] .+ [Δx, Δy] for i in 1:n]))
        @test isempty(_failing(i -> _fwf_padding_is_zero(fwf, i), 1:n))
    end
    # Same as the CPU version, exactly: both add Δx and Δy to every coordinate
    @testset "_translate_floe! on $backend with $FT" for backend in test_backends(),
            FT in (Float64, Float32)
        floes = _make_timestep_test_floes(FT, FloeSettings(FT))
        n = length(floes)
        Δx, Δy = FT(-2e4), FT(1234.5)  # e.g. a ghost through the eastern boundary
        fwf = _launch_on(Subzero._translate_floe!, backend, Subzero.FixedWidthFloes(floes), Δx, Δy)
        foreach(i -> Subzero._translate_floe!(FT, Subzero.get_floe(floes, i), Δx, Δy), 1:n)
        @test isempty(_failing(i -> _fwf_points(fwf, i) == collect(GI.getpoint(floes.poly[i])), 1:n))
        @test isempty(_failing(i -> fwf.centroid[i, :] == floes.centroid[i], 1:n))
        @test isempty(_failing(i -> _fwf_padding_is_zero(fwf, i), 1:n))
    end
end