@testset "Floe" begin
    # Test Setup
    FT = Float64
    Δt = 10
    # test generating monte carlo points
    Lx = 8e4
    Ly = Lx
    Δgrid = 1e4
    hmean = 0.5
    Δh = 0.01
    rng = Xoshiro(1)
    fs_small_min_area = FloeSettings(min_floe_area = 55)

    # Test generating monte carlo points
    file = jldopen("inputs/floe_shapes.jld2", "r")
    floe_coords = file["floe_vertices"][1:end]
    close(file)

    poly1 = Subzero.make_polygon(Subzero.valid_polyvec!(floe_coords[1]))
    centroid1 = GO.centroid(poly1)
    area1 = GO.area(poly1)

    # Test InteractionFields enum
    interactions = range(1, 7)'
    @test interactions[floeidx] == 1
    @test interactions[xpoint] == 4
    @test interactions[overlap] == 7

    # Test with coords inputs
    floe_from_coords = Floe(floe_coords[1], hmean, Δh; u = 0.2, rng = rng)
    @test typeof(floe_from_coords) <: Floe
    @test floe_from_coords.u == 0.2
    @test 0.49 <= floe_from_coords.height <= 0.51
    @test floe_from_coords.centroid == collect(centroid1)
    @test floe_from_coords.area == area1
    @test floe_from_coords.status.tag == Subzero.active
    
    # Test with polygon input
    floe_from_poly = Floe(poly1, hmean, Δh; v = -0.2, rng = Xoshiro(1))
    @test typeof(floe_from_poly) <: Floe
    @test floe_from_poly.u == 0.0
    @test floe_from_poly.v == -0.2
    @test 0.49 <= floe_from_poly.height <= 0.51
    @test floe_from_poly.centroid == collect(centroid1)
    @test floe_from_poly.area == area1
    @test floe_from_poly.status.tag == Subzero.active

    # Test poly_to_floes
    rect_coords = [[[0.0, 0.0], [0.0, 5.0], [10.0, 5.0], [10.0, 0.0], [0.0, 0.0]]]
    rect_poly = Subzero.make_polygon(rect_coords)
    rmax_rect = 2sqrt(5^2 + 2.5^2)
    c_hole_coords = [
        [[0.0, 0.0], [0.0, 10.0], [10.0, 10.0], [10.0, 0.0],[4.0, 0.0], [4.0, 6.0], [2.0, 6.0], [2.0, 0.0], [0.0, 0.0]],
        [[6.0, 4.0], [6.0, 6.0], [7.0, 6.0], [7.0, 4.0], [6.0, 4.0]],
    ]
    c_hole_poly = Subzero.make_polygon(c_hole_coords)
    rmax_cpoly = 2sqrt(5^2 + 5^2)
    # Test polygon with no holes
    floe_arr = StructArray{Floe{FT}}(undef, 0)
    n_new = Subzero.poly_to_floes!(FT, floe_arr, rect_poly, hmean, Δh, rmax_rect)
    @test n_new == 1 && length(floe_arr) == 1
    @test !Subzero.hashole(floe_arr.coords[1])

    # Test with polygon below minimum floe area
    n_new = Subzero.poly_to_floes!(FT, floe_arr, rect_poly, hmean, Δh, rmax_rect; floe_settings = fs_small_min_area)
    @test n_new == 0 && length(floe_arr) == 1

    # Test with polygon with a hole that is split into 3 polyons
    n_new = Subzero.poly_to_floes!(FT, floe_arr, c_hole_poly, hmean, Δh, rmax_cpoly)
    @test n_new == 3 && length(floe_arr) == 4
    @test !any(Subzero.hashole.(floe_arr.coords))

    # Test initialize_floe_field from coord list
    grid = RegRectilinearGrid(; x0 = -Lx, xf = Lx, y0 = -Ly, yf = Ly, Δx = Δgrid, Δy = Δgrid)
    nbound = CollisionBoundary(North; grid)
    sbound = CollisionBoundary(South; grid)
    ebound = CollisionBoundary(East; grid)
    wbound = CollisionBoundary(West; grid)
    domain_no_topo = Domain(; north = nbound, south = sbound, east = ebound, west = wbound)
    island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
    topo1 = [[[-8e4, -8e4], [-8e4, 8e4], [-6e4, 8e4], [-5e4, 4e4], [-6e4, -8e4], [-8e4, -8e4]]]
    domain_with_topo = Domain(;
        north = nbound,
        south = sbound,
        east = ebound,
        west = wbound,
        topography = initialize_topography_field(; coords = [island, topo1])
    )
    topo_polys = Subzero.make_multipolygon([island, topo1])

    # From file without topography -> floes below recommended area
    floe_arr = (@test_logs (:warn, "Some user input floe areas are less than the suggested minimum floe area.") initialize_floe_field(
        FT,
        floe_coords,
        domain_no_topo,
        0.5,
        0.1;
    ))
    nfloes = length(floe_coords)
    @test typeof(floe_arr) <: StructArray{<:Floe}
    @test length(floe_arr) == nfloes
    @test all(floe_arr.id .== range(1, nfloes))

    # From file with small domain -> floes outside of domain
    small_grid = RegRectilinearGrid(; x0 = -Lx/2, xf = Lx/2, y0 = -Ly/2, yf = Ly, Δx = Δgrid, Δy = Δgrid)
    small_domain_no_topo = Domain(;
        north = CollisionBoundary(North; grid = small_grid),
        south = CollisionBoundary(South; grid = small_grid),
        east = CollisionBoundary(East; grid = small_grid),
        west = CollisionBoundary(West; grid = small_grid),
    )
    floe_arr = (@test_logs (:warn, "Some floe centroids are out of the domain.") initialize_floe_field(
        FT,
        floe_coords,
        small_domain_no_topo,
        hmean,
        Δh;
        floe_settings = FloeSettings(min_floe_area = 1e5),
    ))
    @test typeof(floe_arr) <: StructArray{<:Floe}
    @test all(floe_arr.area .> 1e5)

    # From file with topography -> no intersection between topo and floes
    floe_arr = initialize_floe_field(
        FT,
        floe_coords,
        domain_with_topo,
        hmean,
        Δh;
        floe_settings = FloeSettings(min_floe_area = 10),
        rng = Xoshiro(0)
    )
    @test typeof(floe_arr) <: StructArray{<:Floe}
    @test all([sum(GO.area, Subzero.intersect_polys(Subzero.make_polygon(c), topo_polys); init = 0.0) for c in floe_arr.coords] .< 1e-6)
    
    # Test generate_voronoi_coords - general case
    domain_coords = [[[[1, 2], [1.5, 3.5], [1, 5], [2.5, 5], [2.5, 2], [1, 2]]]]
    bounding_box = [[[1, 2], [1, 5], [2.5, 5], [2.5, 2], [1, 2]]]
    voronoi_coords = Subzero.generate_voronoi_coords(
        10,
        [1.5, 3],
        [1, 2],
        domain_coords,
        Xoshiro(1),
        10,
        max_tries = 20, # 20 tries makes it very likely to reach 10 polygons
    )
    bounding_poly = Subzero.make_polygon(bounding_box)
    @test length(voronoi_coords) == 10
    for c in voronoi_coords
        fpoly = Subzero.make_polygon(c)
        @test isapprox(
            sum(GO.area, Subzero.intersect_polys(fpoly, bounding_poly); init = 0.0),
            GO.area(fpoly),
            atol = 1e-3,
        )
    end

    # Test warning and no points generated
    warning_str = "Only 0 floes were able to be generated in 10 tries during \
        voronoi tesselation."
    @test @test_logs (:warn, warning_str) Subzero.generate_voronoi_coords(
        0, # Don't generate any points
        [1.5, 3],
        [1, 2],
        domain_coords,
        Xoshiro(1),
        10, # Higher min to warn so we can test warning
    ) == Vector{Vector{Vector{Float64}}}()
    
    # Test initialize_floe_field with voronoi
    floe_arr = initialize_floe_field(
        FT,
        25,
        [0.5],
        domain_with_topo,
        0.5,
        0.1;
        floe_settings = FloeSettings(min_floe_area = 1e4),
        rng = Xoshiro(1)
    )
    @test isapprox(
        sum(floe_arr.area)/(1.6e5*1.6e5 - GO.area(topo_polys)),
        0.5,
        atol = 1e-1
    )
    @test all(floe_arr.area .> 1e4)
    for (i, c) in enumerate(floe_arr.coords)
        Subzero.intersect_polys(Subzero.make_polygon(c), topo_polys)
    end
    @test all([sum(GO.area, Subzero.intersect_polys(Subzero.make_polygon(c), topo_polys); init = 0.0) for c in floe_arr.coords] .< 1e-6)

    nfloes = length(floe_arr)
    @test all(floe_arr.id .== range(1, nfloes))

    concentrations = [1 0.3; 0 0.5]
    rng = Xoshiro(2)
    floe_arr = initialize_floe_field(
        25,
        concentrations,
        domain_with_topo,
        0.5,
        0.1;
        floe_settings = FloeSettings(min_floe_area = 1e4),
        rng = rng
    )
    nfloes = length(floe_arr)
    floe_polys = [Subzero.make_polygon(f) for f in floe_arr.coords]
    first_cell = [[[-8e4, -8e4], [-8e4, 0], [0, 0], [0, -8e4], [-8e4, -8e4]]]
    for j in 1:2
        for i in 1:2
            cell = Subzero.make_polygon(Subzero.translate(first_cell, 8e4*(j-1), 8e4*(i-1)))
            cell_without_topos = Subzero.diff_polys(cell, topo_polys)
            open_cell_area = sum(GO.area, cell_without_topos; init = 0.0)
            c = concentrations[i, j]
            floes_in_cell_area = 0
            n_polys = 0
            for floe in floe_polys
                for mask in cell_without_topos
                    floes_in_cell_area += GO.area(Subzero.intersect_polys(floe, mask))
                end
            end
            @test c - 100eps() <= floes_in_cell_area/open_cell_area < 1 + eps()
        end
    end
    @test all([sum(GO.area, Subzero.intersect_polys(p, topo_polys); init = 0.0) for p in floe_polys] .< 1e-3)
    @test all(floe_arr.id .== range(1, nfloes))
    @test typeof(initialize_floe_field(
        Float32,
        25,
        concentrations,
        domain_with_topo,
        0.5,
        0.1;
    )) <: StructArray{<:Floe{Float32}}
end