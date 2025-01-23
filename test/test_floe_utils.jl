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

    # Test predicate hashole for coords, polygons and multipolygons
    @test !Subzero.hashole([ext])
    @test Subzero.hashole([ext, hole1])
    @test !Subzero.hashole(poly_nohole)
    @test Subzero.hashole(poly_hole1)
    @test Subzero.hashole(poly_hole2)

    # Test removing holes from polygons
    copy_holes = [ext, hole1]
    poly_copy_holes = Subzero.make_polygon(copy_holes)
    Subzero.rmholes!(copy_holes)
    @test copy_holes == [ext]
    Subzero.rmholes!(poly_copy_holes)
    @test GI.nhole(poly_copy_holes) == 0

    # Test translating coordinates and polygons
    @test Subzero.translate([ext], 0.0, 0.0) == [ext]
    trans_ext = Subzero.translate([ext], 1.0, 2.0)
    @test trans_ext == [[[1.0, 3.0],  [1.0, 2.0],  [2.0, 2.0],
                        [2.0, 3.0], [1.0, 3.0]]]
    copy_ext = [deepcopy(ext)]
    Subzero.translate!(copy_ext, 1.0, 2.0)
    @test copy_ext == trans_ext
    test_trans = [[[-2.0, 2.0], [-2.0, 1.0], [-1.0, 1.0], [-1.0, 2.0]]]
    @test Subzero.translate(test_trans, 1.5, -1.5) ==
        [[[-0.5, 0.5], [-0.5, -0.5], [0.5, -0.5], [0.5, 0.5]]]
    Subzero.translate!(test_trans, 1.5, -1.5)
    @test test_trans == [[[-0.5, 0.5], [-0.5, -0.5], [0.5, -0.5], [0.5, 0.5]]]

    # Test moment of intertia calculations - compared to values output my MATLAB
    poly_moment = Subzero._calc_moment_inertia(Float64, Subzero.make_polygon([ext]), [0.5, 0.5], 0.25)
    @test isapprox(poly_moment, 38.333, atol = 0.001)
    @test Subzero._calc_moment_inertia(Float64, poly_nohole, GO.centroid(poly_nohole), 0.25) == poly_moment
    tri_poly = Subzero.make_polygon([[[0, 1], [0, 0], [1, 0], [0, 1]]] .* 6.67)
    tri_moment = Subzero._calc_moment_inertia(Float64, tri_poly, GO.centroid(tri_poly), 0.5)
    @test isapprox(tri_moment, 50581.145, atol = 0.001)

    # ------------------------- Test finding shared points ---------------------
    two_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [20.0, 0.0], [0.0, 0.0]]],
        [[[20.0, 0.0], [20.0, 20.0], [40.0, 20.0], [40.0, 0.0], [20.0, 0.0]]]
    ]
    @test Subzero.which_vertices_match_points(
        two_shared_v[1][1],
        Subzero.make_polygon(two_shared_v[2]),
    ) == [1, 2]

    @test Subzero.which_points_on_edges(
        two_shared_v[1][1],
        two_shared_v[2],
    ) == [3, 4]

    three_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [20.0, 10.0], [20.0, 0.0], [0.0, 0.0]]],
        [[[40.0, 20.0], [40.0, 0.0], [20.0, 0.0], [20.0, 10.0], [20.0, 20.0], [40.0, 20.0]]]
    ]
    @test Subzero.which_vertices_match_points(
        three_shared_v[1][1],
        Subzero.make_polygon(three_shared_v[2]),
    ) == [3, 4, 5]

    @test Subzero.which_points_on_edges(
        three_shared_v[1][1],
        three_shared_v[2],
    ) == [3, 4, 5]

    deleteat!(three_shared_v[2][1], 4)
    @test Subzero.which_points_on_edges(
        three_shared_v[1][1],
        three_shared_v[2],
    ) == [3, 4, 5]

    four_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [20.0, 18.0], [20.0, 15.0], [20.0, 0.0], [0.0, 0.0]]],
        [[[20.0, 18.0], [20.0, 20.0], [40.0, 20.0], [40.0, 0.0], [20.0, 0.0], [20.0, 15.0], [20.0, 18.0]]]
    ]
    @test Subzero.which_vertices_match_points(
        four_shared_v[1][1],
        Subzero.make_polygon(four_shared_v[2]),
    ) == [1, 2, 5, 6]

    offset_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [20.0, 0.0], [0.0, 0.0]]],
        [[[5.0, 20.0], [5.0, 25.0], [25.0, 25.0], [25.0, 20.0], [10.0, 20.0], [5.0, 20.0]]]
    ]
    @test Subzero.which_points_on_edges(
        offset_shared_v[1][1],
        offset_shared_v[2],
    ) == [3]
    @test Subzero.which_points_on_edges(
        offset_shared_v[2][1],
        offset_shared_v[1],
    ) == [1, 5]

    triange_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [5.0, 5.0], [0.0, 0.0]]],
        [[[0.0, 0.0], [5.0, 5.0], [20.0, 20.0], [20.0, 0.0], [0.0, 0.0]]]
    ]
    @test Subzero.which_vertices_match_points(
        triange_shared_v[1][1],
        Subzero.make_polygon(triange_shared_v[2]),
    ) == [1, 2, 3]

    @test Subzero.which_points_on_edges(
        triange_shared_v[1][1],
        triange_shared_v[2],
    ) == [1, 3, 4]

    deleteat!(triange_shared_v[2][1], 2)
    @test Subzero.which_points_on_edges(
        triange_shared_v[1][1],
        triange_shared_v[2],
    ) == [1, 3, 4]

    squares_midpoint = (20.0, 10.0)
    triangle_midpoint = (10.0, 10.0)
    @test Subzero.find_shared_edges_midpoint(
        two_shared_v[1],
        two_shared_v[2],
    ) == squares_midpoint
    @test Subzero.find_shared_edges_midpoint(
        three_shared_v[1],
        three_shared_v[2],
    ) == squares_midpoint
    @test Subzero.find_shared_edges_midpoint(
        four_shared_v[1],
        four_shared_v[2],
    ) == squares_midpoint
    @test Subzero.find_shared_edges_midpoint(
        triange_shared_v[1],
        triange_shared_v[2],
    ) == (10.0, 10.0)
    @test Subzero.find_shared_edges_midpoint(
        offset_shared_v[2],
        offset_shared_v[1],
    ) == (7.5, 20)

    # ------------------ Test rotating coordinates ------------------
    og_coords = [
        [[-1.0, -1.0], [-1, 1], [1, 1], [1, -1], [-1, -1]],
        [[-1.0, -1.0], [-1, 1], [1, 1], [1, -1], [-1, -1]],
    ]
    copy_coords = deepcopy(og_coords)
    Subzero.rotate_radians!(copy_coords, π/4)
    same = true
    answer = [[[0, -√2], [-√2, 0], [0, √2], [√2, 0], [0, -√2]],
        [[0, -√2], [-√2, 0], [0, √2], [√2, 0], [0, -√2]]]
    for i in eachindex(copy_coords)
        same = same && all(isapprox.(copy_coords[i], answer[i]))
    end
    @test same
    Subzero.rotate_radians!(copy_coords, 7π/4)
    same = true
    for i in eachindex(copy_coords)
        same = same & all(isapprox.(copy_coords[i], og_coords[i]))
    end
    @test same
end