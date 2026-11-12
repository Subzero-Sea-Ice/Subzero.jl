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
    fs_no_min_area = FloeSettings(min_floe_area = 0)

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
    floe_from_coords = Floe(floe_coords[1], hmean; u = 0.2, rng = rng)
    @test typeof(floe_from_coords) <: Floe
    @test floe_from_coords.u == 0.2
    @test floe_from_coords.height == 0.5
    @test floe_from_coords.centroid == collect(centroid1)
    @test floe_from_coords.area == area1
    @test floe_from_coords.status.tag == Subzero.active
    
    # Test with polygon input
    floe_from_poly = Floe(poly1, hmean; v = -0.2, rng = Xoshiro(1))
    @test typeof(floe_from_poly) <: Floe
    @test floe_from_poly.u == 0.0
    @test floe_from_poly.v == -0.2
    @test floe_from_poly.height == 0.5
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
    n_new = Subzero._poly_to_floes!(FT, floe_arr, rect_poly, hmean, Δh, rmax_rect; floe_settings = fs_no_min_area)
    @test n_new == 1 && length(floe_arr) == 1
    @test !Subzero.hashole(floe_arr.poly[1])

    # Test with polygon below minimum floe area
    n_new = Subzero._poly_to_floes!(FT, floe_arr, rect_poly, hmean, Δh, rmax_rect; floe_settings = fs_small_min_area)
    @test n_new == 0 && length(floe_arr) == 1

    # Test with polygon with a hole that is split into 3 polyons
    n_new = Subzero._poly_to_floes!(FT, floe_arr, c_hole_poly, hmean, Δh, rmax_cpoly; floe_settings = fs_no_min_area)
    @test n_new == 3 && length(floe_arr) == 4
    @test !any(Subzero.hashole.(floe_arr.poly))

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
        floe_settings = fs_no_min_area,
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
    @test all([sum(GO.area, Subzero.intersect_polys(p, topo_polys); init = 0.0) for p in floe_arr.poly] .< 1e-6)
    
    # Test _generate_voronoi_coords - general case
    domain_coords = [[[[1, 2], [1.5, 3.5], [1, 5], [2.5, 5], [2.5, 2], [1, 2]]]]
    domain_poly = Subzero.make_multipolygon(domain_coords)
    bounding_box = [[[1, 2], [1, 5], [2.5, 5], [2.5, 2], [1, 2]]]
    bounding_poly = Subzero.make_polygon(bounding_box)
    voronoi_polys = Subzero._generate_voronoi_coords(
        FT,
        10,
        1.5, 3, # Δx, Δy
        1, 2, # xmin, ymin
        domain_poly,
        Xoshiro(1),
        10,
        max_tries = 20, # 20 tries makes it very likely to reach 10 polygons
    )
    @test length(voronoi_polys) == 10
    for fpoly in voronoi_polys
        @test isapprox(
            sum(GO.area, Subzero.intersect_polys(fpoly, bounding_poly); init = 0.0),
            GO.area(fpoly),
            atol = 1e-3,
        )
    end

    # Test warning and no points generated
    warning_str = "Only 0 floes were able to be generated in 10 tries during \
        voronoi tesselation."
    @test @test_logs (:warn, warning_str) Subzero._generate_voronoi_coords(FT,
        0, # Don't generate any points
        1.5, 3, # Δx, Δy
        1, 2, # xmin, ymin
        domain_poly,
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
    for (i, p) in enumerate(floe_arr.poly)
        Subzero.intersect_polys(p, topo_polys)
    end
    @test all([sum(GO.area, Subzero.intersect_polys(p, topo_polys); init = 0.0) for p in floe_arr.poly] .< 1e-6)

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
    floe_polys = floe_arr.poly
    first_cell = [[[-8e4, -8e4], [-8e4, 0], [0, 0], [0, -8e4], [-8e4, -8e4]]]
    for j in 1:2
        for i in 1:2
            cell = Subzero.make_polygon(translate_coords(first_cell, 8e4*(j-1), 8e4*(i-1)))
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
        floe_settings = fs_no_min_area,
    )) <: StructArray{<:Floe{Float32}}

    @testset "FixedWidthFloes" begin
        FT = Float64
        floes = _make_timestep_test_floes(FT, FloeSettings())
        fwf = Subzero.FixedWidthFloes(floes)
        n = length(floes)
        # Ragged fields are padded to the longest floe
        n_points = [GI.npoint(p) for p in floes.poly]
        @test fwf.n_points[1:n] == n_points
        @test size(fwf.poly) == (Subzero.capacity(fwf), maximum(n_points), 2)
        @test fwf.num_inters[1:n] == floes.num_inters
        @test size(fwf.interactions) == (Subzero.capacity(fwf), maximum(floes.num_inters), 7)
        for i in eachindex(floes)
            @test all(fwf.poly[i, (n_points[i] + 1):end, :] .== 0)
            @test [(fwf.poly[i, j, 1], fwf.poly[i, j, 2]) for j in 1:n_points[i]] ==
                collect(GI.getpoint(floes.poly[i]))
            ninters = floes.num_inters[i]
            @test fwf.interactions[i, 1:ninters, :] == floes.interactions[i][1:ninters, :]
            @test all(fwf.interactions[i, (ninters + 1):end, :] .== 0)
            @test fwf.collision_force[i, :] == vec(floes.collision_force[i])
            @test fwf.stress_accum[i, :, :] == floes.stress_accum[i]
        end
        @test fwf.rmax[1:n] == floes.rmax
        @test fwf.status[1:n] == fill(Int32(Subzero.active), n)
        @test all(iszero, fwf.ghost_id) && all(iszero, fwf.num_ghosts) && all(iszero, fwf.ghosts)
        # Ghosts, ghost ids, and status
        ghosted = deepcopy(floes)
        ghosted.ghosts[1] = [8, 9, 10]
        ghosted.ghosts[2] = [7]
        ghosted.ghost_id[7:10] .= [1, 1, 2, 3]
        ghosted.status[3].tag = Subzero.remove
        ghosted.status[4].tag = Subzero.fuse
        fwf = Subzero.FixedWidthFloes(ghosted)
        @test size(fwf.ghosts) == (Subzero.capacity(fwf), Subzero.MAX_GHOSTS)
        @test fwf.num_ghosts[1:n] == length.(ghosted.ghosts)
        for i in eachindex(ghosted)
            @test fwf.ghosts[i, 1:fwf.num_ghosts[i]] == ghosted.ghosts[i]
            @test all(iszero, fwf.ghosts[i, (fwf.num_ghosts[i] + 1):end])
        end
        @test fwf.ghost_id[1:n] == ghosted.ghost_id
        @test fwf.status[1:4] == Int32.([Subzero.active, Subzero.active, Subzero.remove, Subzero.fuse])
        # More ghosts than MAX_GHOSTS still fit
        ghosted.ghosts[1] = [7, 8, 9, 10]
        @test Subzero.FixedWidthFloes(ghosted).ghosts[1, :] == [7, 8, 9, 10]
        # Does not share memory with floes
        fwf.height .= 0
        fwf.centroid .= 0
        @test all(floes.height .> 0)
        @test all(c -> all(c .> 0), floes.centroid)
        # Round trip leaves floes unchanged
        floes = _make_timestep_test_floes(FT, FloeSettings())
        expected = deepcopy(floes)
        Subzero.update_floes!(floes, Subzero.adapt(Array, Subzero.FixedWidthFloes(floes)))
        for field in propertynames(floes)
            field == :status && continue  # Status has no ==, compared below
            @test getproperty(floes, field) == getproperty(expected, field)
        end
        @test [s.tag for s in floes.status] == [s.tag for s in expected.status]
        # Changes are written back to floes
        fwf = Subzero.FixedWidthFloes(floes)
        fwf.poly[3, 1:fwf.n_points[3], :] .+= 1
        fwf.centroid[3, :] .+= 1
        fwf.strain[3, 1, 2] = 42
        fwf.u[3] = 42
        Subzero.update_floes!(floes, fwf)
        @test collect(GI.getpoint(floes.poly[3])) ==
            [(x + 1, y + 1) for (x, y) in GI.getpoint(expected.poly[3])]
        @test floes.centroid[3] == expected.centroid[3] .+ 1
        @test floes.strain[3][1, 2] == 42
        @test floes.u[3] == 42
        @test floes.poly[3] isa Subzero.Polys{FT}
    end

    @testset "FixedWidthFloes with spare capacity" begin
        FT = Float64
        floes = _make_timestep_test_floes(FT, FloeSettings())
        n = length(floes)
        # By default 25% more room than floes, rounded up
        @test Subzero.FixedWidthFloes(floes).n_floes == n
        @test Subzero.capacity(Subzero.FixedWidthFloes(floes)) == ceil(Int, 1.25n)
        @test Subzero.capacity(Subzero.FixedWidthFloes(floes[1:0])) == 0
        @test Subzero.capacity(Subzero.FixedWidthFloes(floes[1:1])) == 2
        fwf = Subzero.FixedWidthFloes(floes; capacity = n)
        @test fwf.n_floes == Subzero.capacity(fwf) == n
        spare = Subzero.FixedWidthFloes(floes; capacity = n + 3)
        @test spare.n_floes == n
        @test Subzero.capacity(spare) == n + 3
        # The floes in use are the same as without spare capacity, the spare rows are zero
        for field in fieldnames(Subzero.FixedWidthFloes)
            field == :n_floes && continue
            a, b = getfield(spare, field), getfield(fwf, field)
            @test size(a) == (n + 3, size(b)[2:end]...)
            @test selectdim(a, 1, 1:n) == b
            @test all(iszero, selectdim(a, 1, (n + 1):(n + 3)))
        end
        @test_throws ArgumentError Subzero.FixedWidthFloes(floes; capacity = n - 1)
        # Round trip with spare capacity leaves floes unchanged
        expected = deepcopy(floes)
        Subzero.update_floes!(floes, spare)
        @test floes.centroid == expected.centroid
        @test floes.u == expected.u
        @test floes.poly == expected.poly
        # Changing the number of floes in use shares the arrays
        fewer = Subzero.set_n_floes(spare, n - 2)
        @test fewer.n_floes == n - 2
        @test Subzero.capacity(fewer) == n + 3
        @test fewer.poly === spare.poly && fewer.flags === spare.flags
        @test Subzero.set_n_floes(spare, 0).n_floes == 0
        @test Subzero.set_n_floes(spare, n + 3).n_floes == n + 3
        @test_throws ArgumentError Subzero.set_n_floes(spare, n + 4)
        @test_throws ArgumentError Subzero.set_n_floes(spare, -1)
        @test_throws ArgumentError Subzero.update_floes!(floes, fewer)
        #= Kernels only run on the floes in use: mark! writes `floe_ran` into the flags of
        every floe it runs on, and the other floes keep the `floe_did_not_run` flags they
        start with. These are test markers, not FLAG_* event bits. =#
        floe_ran = 0x01
        floe_did_not_run = 0x00
        mark!(floes, i, value) = (floes.flags[i] = value; nothing)
        @testset "launch_per_floe! on $backend" for backend in test_backends()
            for m in (0, n - 2, n + 3)
                # A new copy each time: adapt to the CPU backend shares the arrays
                fresh = Subzero.FixedWidthFloes(floes; capacity = n + 3)
                dev_floes = Subzero.adapt(backend, Subzero.set_n_floes(fresh, m))
                @test dev_floes.n_floes == m
                Subzero.launch_per_floe!(mark!, backend, dev_floes, floe_ran)
                Subzero.KernelAbstractions.synchronize(backend)
                host_floes = Subzero.adapt(Array, dev_floes)
                @test host_floes.n_floes == m
                @test host_floes.flags ==
                    [i <= m ? floe_ran : floe_did_not_run for i in 1:(n + 3)]
            end
        end
    end
end
