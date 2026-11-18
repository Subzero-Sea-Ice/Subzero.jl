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
end