@testset "Coupling" begin
    FT = Float64
    @testset "Generating sub-floe points" begin
        FT = Float64
        # test generating monte carlo points
        file = jldopen("inputs/floe_shapes.jld2", "r")
        floe_coords = file["floe_vertices"][1:end]
        close(file)
        poly1 = Subzero.make_polygon(Subzero.valid_polyvec!(floe_coords[1]))
        centroid1 = GO.centroid(poly1)
        origin_coords = Subzero.translate(
            floe_coords[1],
            -centroid1[1],
            -centroid1[2],
        )
        origin_poly = Subzero.make_polygon(origin_coords)
        origin_centroid = GO.centroid(origin_poly)
        xo, yo = first.(origin_coords[1]), last.(origin_coords[1])
        rmax = sqrt(maximum([sum(xo[i]^2 + yo[i]^2) for i in eachindex(xo)]))
        area = GO.area(poly1)
        mc_x, mc_y, status = Subzero.generate_subfloe_points(
            MonteCarloPointsGenerator(),
            origin_poly,
            origin_centroid,
            area,
            Subzero.Status(),
            Xoshiro(1)
        )
        @test length(mc_x) == length(mc_y) && length(mc_x) > 0
        mc_in = [GO.coveredby((mc_x[i], mc_y[i]), origin_poly) for i in eachindex(mc_x)]
        @test all(mc_in)
        xmin, xmax = extrema(xo)
        ymin, ymax = extrema(yo)
        @test abs(sum(mc_in)/1000 * (xmax - xmin) * (ymax - ymin) - area)/area < 0.1
        @test status.tag == Subzero.active
        # Test that random number generator is working
        mc_x2, mc_y2, status2 = Subzero.generate_subfloe_points(
            MonteCarloPointsGenerator(),
            origin_poly,
            origin_centroid,
            area,
            Subzero.Status(),
            Xoshiro(1)
        )
        @test all(mc_x .== mc_x2)
        @test all(mc_y .== mc_y2)
        @test status2.tag == Subzero.active
    
        mc_x3, mc_y3, status3 = Subzero.generate_subfloe_points(
            MonteCarloPointsGenerator{Float32}(),
            GO.tuples(origin_poly, Float32),
            origin_centroid,
            area,
            Subzero.Status(),
            Xoshiro(1)
        )
        @test status3.tag == Subzero.active
        @test eltype(mc_x3) == eltype(mc_y3) == Float32

        # test generating sub-grid points for grid with Δx = Δy = 10
        point_generator = SubGridPointsGenerator{Float64}(10/sqrt(2))
        # Floe is smaller than grid cells --> centroid and vertices added
        square = Subzero.make_polygon([[
            [-2.5, -2.5],
            [-2.5, 2.5],
            [2.5, 2.5],
            [2.5, -2.5],
            [-2.5, -2.5],
        ]])
        square_centroid = GO.centroid(square)
        xpoints, ypoints = Subzero.generate_subfloe_points(
            point_generator,
            square,
            square_centroid,
            0.0, # Not used
            Subzero.Status(),
            Xoshiro(), # Not random
        )
        @test xpoints == [-2.5, -2.5, 2.5, 2.5, 0.0]
        @test ypoints == [-2.5, 2.5, 2.5, -2.5, 0.0]
        # Floe is larger than grid cell
        tall_rect = Subzero.make_polygon([[
            [-2.0, -10.0],
            [-2.0, 10.0],
            [2.0, 10.0],
            [2.0, -10.0],
            [-2.0, -10.0],
        ]])
        tall_rect_centroid = GO.centroid(tall_rect)
        xpoints, ypoints = Subzero.generate_subfloe_points(
            point_generator,
            tall_rect,
            tall_rect_centroid,
            0.0, # Not used
            Subzero.Status(),
            Xoshiro(), # Not random
        )
        @test xpoints == [repeat([-2.0], 5); repeat([2.0], 5); repeat([0.0], 3)]
        @test all(isapprox.(
            ypoints,
            [-10.0, -6.46447, 0.0, 6.46447, 10.0, 10.0, 6.46447, 0.0,
                -6.46447, -10, -6.46447, 0.0, 6.46447,
            ],
            atol = 1e-5))

        wide_rect = Subzero.make_polygon([[
            [-10.0, -2.0],
            [-10.0, 2.0],
            [10.0, 2.0],
            [10.0, -2.0],
            [-10.0, -2.0],
        ]])
        wide_rect_centroid = GO.centroid(wide_rect)
        xpoints, ypoints = Subzero.generate_subfloe_points(
            point_generator,
            wide_rect,
            wide_rect_centroid,
            0.0, # Not used
            Subzero.Status(),
            Xoshiro(), # Not random
        )
        @test all(isapprox.(
            xpoints,
            [-10, -10, -6.46447, 0.0, 6.46447, 10, 10, 6.46447, 0.0,
                -6.464466, -6.46447, 0, 6.46447,
            ],
            atol = 1e-5
        ))
        @test ypoints == [-2; repeat([2], 5); repeat([-2], 4); repeat([0], 3)] 
        trapeziod = Subzero.make_polygon([[
            [-8.0, -8.0],
            [-4.0, 8.0],
            [4.0, 8.0],
            [8.0, -8.0],
            [-8.0, -8.0],
        ]])
        trapeziod_centroid = GO.centroid(trapeziod)
        xpoints, ypoints = Subzero.generate_subfloe_points(
            point_generator,
            trapeziod,
            trapeziod_centroid,
            0.0, # Not used
            Subzero.Status(),
            Xoshiro(), # Not random
        )
        ypoints .+= trapeziod_centroid[2]  # y-point of centroid is not centered on the origin
        @test all(isapprox.(
            xpoints,
            [-8, -7.14251, -6.0, -4.85749, -4.0, 0.0, 4.0, 4.85749, 6.0,
                7.14251, 8.0, 4.46447, 0.0, -4.46447, -4.46447, 0.0, 4.46447,
                -4.46447, 0.0, 4.46447, -4.46447, 0.0, 4.46447
            ],
            atol = 1e-5
        ))
        @test all(isapprox.(
            ypoints,
            [-8; -4.57003; 0.0; 4.57003; repeat([8.0], 3); 4.57003; 0.0;
                -4.57003; repeat([-8.0], 4); repeat([-4.46447], 3);
                repeat([0.0], 3); repeat([4.46447], 3);
            ],
            atol = 1e-5
        ))
    end

    @testset "Coupling Helper Functions" begin
        grid = Subzero.RegRectilinearGrid(; x0 = -10, xf = 10, y0 = -8, yf = 8, Δx = 2, Δy = 4)
        # Test find_cell_indices
        xpoints = [-10.5, -10, -10, -6.5, -6, -4, 10, 10.5, 12]
        ypoints = [0.0, 6.0, -8.0, 4.5, 0.0, 5.0, -8.0, 0.0, 0.0]
        xidx = [1, 1, 1, 3, 3, 4, 11, 11, 12]
        yidx = [3, 5, 1, 4, 3, 4, 1, 3, 3]
        for i in eachindex(xpoints)
            @test  Subzero.find_center_cell_index(
                xpoints[i],
                ypoints[i],
                grid,
            ) == (xidx[i], yidx[i])
        end

        # Test filter_oob_points
        open_bound = Subzero.OpenBoundary(East; grid)
        periodic_bound = Subzero.PeriodicBoundary(East; grid)
        x = [-12, -10, -8, -6, 0, 4, 4, 10, 12, 12]
        y = [5, -6, 4, 10, -10, 8, -8, -6, 4, 10]
        open_open_answers = [false, true, true, false, false, true, true, true, false, false]
        periodic_open_answers = [false, true, true, true, true, true, true, true, false, false]
        open_periodic_answers = [true, true, true, false, false, true, true, true, true, false]
        periodic_periodic_answers = [true, true, true, true, true, true, true, true, true, true]
        # both bounds non-periodic
        for i in eachindex(x)
            @test Subzero.in_bounds(x[i], y[i], grid, open_bound, open_bound) == open_open_answers[i]
            @test Subzero.in_bounds(x[i], y[i], grid, open_bound, periodic_bound) == open_periodic_answers[i]
            @test Subzero.in_bounds(x[i], y[i], grid, periodic_bound, open_bound) == periodic_open_answers[i]
            @test Subzero.in_bounds(x[i], y[i], grid, periodic_bound, periodic_bound) == periodic_periodic_answers[i]
        end

        # Test find_interp_knots
        # in bounds
        @test Subzero.find_interp_knots(
            [4],
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            2,
            open_bound,
        ) == (0:10:60, 1:7)
        @test Subzero.find_interp_knots(
            [4],
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            2,
            periodic_bound,
        ) == (0:10:60, 1:7)
        # out of bounds to left
        @test Subzero.find_interp_knots(
            [0, 1],
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            2,
            open_bound,
        ) == (0:10:30, 1:4)
        @test Subzero.find_interp_knots(
            [0, 1],
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            2,
            periodic_bound,
        ) == (-40:10:30, [5:8; 1:4])
        # out of bounds to right
        @test Subzero.find_interp_knots(
            [8, 9],
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            1,
            open_bound,
        ) == (50:10:80, 6:9)
        @test Subzero.find_interp_knots(
            [8, 9],
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            1,
            periodic_bound,
        ) == (50:10:100, [6:8; 1:3])
        # buffer out of bounds on both sides
        @test Subzero.find_interp_knots(
            1:8,
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            2,
            open_bound,
        ) == (0:10:80, 1:9)
        @test Subzero.find_interp_knots(
            1:8,
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            2,
            periodic_bound,
        ) == (-30:10:100, [6:8; 1:8; 1:3])
        # floe out of bounds on both sides - would be filtered out for periodic in advance
        @test Subzero.find_interp_knots(
            0:9,
            8,  # 8 grid cells
            0:10:80,  # grid lines
            80, # grid dimension is 80m
            2,
            open_bound,
        ) == (0:10:80, 1:9)

        #Test cell_coords
         @test GO.equals(Subzero.center_cell_coords(Float64, 2, 3, grid, periodic_bound, periodic_bound),
            GI.Polygon([[[-9, -2], [-9, 2], [-7, 2], [-7, -2], [-9, -2]]]))
        @test GO.equals(Subzero.center_cell_coords(Float64, 1, 1, grid, open_bound, open_bound),
            GI.Polygon([[[-10, -8], [-10, -6], [-9, -6], [-9, -8], [-10, -8]]]))
        @test GO.equals(Subzero.center_cell_coords(Float64, 11, 6, grid, periodic_bound, periodic_bound),
            GI.Polygon([[[9, 10], [9, 14], [11, 14], [11, 10], [9, 10]]]))
        @test GO.equals(Subzero.center_cell_coords(Float64, 11, 6, grid, open_bound, open_bound),
            GI.Polygon([[[9, 8], [9, 8], [10, 8], [10, 8], [9, 8]]]))
        @test GO.equals(Subzero.center_cell_coords(Float64, 11, 6, grid, open_bound, periodic_bound),
            GI.Polygon([[[9, 8], [9, 8], [11, 8], [11, 8], [9, 8]]]))
        @test GO.equals(Subzero.center_cell_coords(Float64, 11, 6, grid, periodic_bound, open_bound),
            GI.Polygon([[[9, 10], [9, 14], [10, 14], [10, 10], [9, 10]]]))

        # Test aggregate_grid_stress!
        function test_floe_to_grid(
            floeidx,
            xidx,
            yidx,
            τxs,
            τys,
            grid,
            bound1,
            bound2,
            scells,
            settings,
            occupied_cells,
            Δxs,
            Δys,
            sum_τx,
            sum_τy,
            npoints,
        )
            for i in eachindex(xidx)
                Subzero.floe_to_grid_info!(
                    floeidx,
                    xidx[i],
                    yidx[i],
                    τxs[i],
                    τys[i],
                    grid,
                    bound1,
                    bound2,
                    scells,
                    settings,
                )
            end

            for i in eachindex(occupied_cells)
                @test grid.floe_locations[occupied_cells[i]].floeidx[end] == floeidx
                @test grid.floe_locations[occupied_cells[i]].Δx[end] == Δxs[i]
                @test grid.floe_locations[occupied_cells[i]].Δy[end] == Δys[i]
                @test scells[occupied_cells[i]].τx[end] == sum_τx[i]
                @test scells[occupied_cells[i]].τy[end] == sum_τy[i]
                @test scells[occupied_cells[i]].npoints[end] == npoints[i]
            end

            empty_check = true
            for cidx in CartesianIndices(scells)
                if !(cidx in occupied_cells)
                    empty_check = empty_check &&
                        isempty(grid.floe_locations[cidx].floeidx) &&
                        isempty(grid.floe_locations[cidx].Δx) &&
                        isempty(grid.floe_locations[cidx].Δy)
                    empty_check = empty_check &&
                        isempty(scells[cidx].τx) &&
                        isempty(scells[cidx].τy) &&
                        isempty(scells[cidx].npoints)
                end
            end
            @test empty_check
        end
        # Open bounds, everything within bounds
        test_floe_to_grid(
            1,
            [7, 7, 6, 6, 7, 7],
            [4, 4, 3, 3, 4, 4],
            ones(Float64, 6),
            2ones(Float64, 6),
            deepcopy(grid),
            open_bound,
            open_bound,
            Subzero.Ocean(; grid, u = 0, v = 0, temp = 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [CartesianIndex(7, 4), CartesianIndex(6, 3)],
            [0.0, 0.0],
            [0.0, 0.0],
            [-4, -2],
            [-8, -4],
            [4, 2],
        )
        # Periodic bounds, everything within bounds
        test_floe_to_grid(
            2,
            [7, 7, 8, 8, 9, 9],
            [2, 3, 3, 3, 2, 2],
            ones(Float64, 6),
            2ones(Float64, 6),
            deepcopy(grid),
            periodic_bound,
            periodic_bound,
            Subzero.Ocean(; grid, u = 0, v = 0, temp = 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [
                CartesianIndex(7, 2),
                CartesianIndex(7, 3),
                CartesianIndex(9, 2),
                CartesianIndex(8, 3),
            ],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [-1, -1, -2, -2],
            [-2, -2, -4, -4],
            [1, 1, 2, 2],
        )
        # Periodic bounds in NS, floe out of bounds - moved to other side of grid
        test_floe_to_grid(
            3,
            [10, 10, 10, 11, 11, 11, 11],
            [4, 5, 6, 5, 6, 5, 6],
            ones(Float64, 7),
            2ones(Float64, 7),
            deepcopy(grid),
            periodic_bound,
            open_bound,
            Subzero.Ocean(; grid, u = 0, v = 0, temp = 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [
                CartesianIndex(10, 1),
                CartesianIndex(11, 1),
                CartesianIndex(10, 2),
                CartesianIndex(11, 2),
                CartesianIndex(10, 4),
            ],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [-16.0, -16.0, -16.0, -16.0, 0.0],
            [-1, -2, -1, -2, -1],
            [-2, -4, -2, -4, -2],
            [1, 2, 1, 2, 1],
        )
        # Periodic bounds in EW, floe out of bounds - moved to other side of grid
        test_floe_to_grid(
            4,
            [11, 11, 12, 12, 11],
            [4, 5, 5, 5, 4],
            ones(Float64, 5),
            2ones(Float64, 5),
            deepcopy(grid),
            open_bound,
            periodic_bound,
            Subzero.Ocean(; grid, u = 0, v = 0, temp = 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [
                CartesianIndex(1, 4),
                CartesianIndex(1, 5),
                CartesianIndex(2, 5),
            ],
            [-20.0, -20.0, -20.0],
            [0.0, 0.0, 0.0],
            [-2, -1, -2],
            [-4, -2, -4],
            [2, 1, 2],
        )
        # Periodic bounds in both direction, floe out of bounds
        test_floe_to_grid(
            2,
            [0, -1, -1, 1, -1],
            [0, -1, -2, 1, -1],
            -ones(Float64, 5),
            -2ones(Float64, 5),
            deepcopy(grid),
            periodic_bound,
            periodic_bound,
            Subzero.Ocean(; grid, u = 0, v = 0, temp = 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [
                CartesianIndex(1, 1),
                CartesianIndex(10, 4),
                CartesianIndex(9, 3),
                CartesianIndex(9, 2),
            ],
            [0.0, 20.0, 20.0, 20.0],
            [0.0, 16.0, 16.0, 16.0],
            [1, 1, 2, 1],
            [2, 2, 4, 2],
            [1, 1, 2, 1],
        )
    end
    @testset "OA Forcings" begin
    #     # set up model and floe
        FT = Float64
        grid = Subzero.RegRectilinearGrid(; x0 = -1e5, xf = 1e5, y0 = -1e5, yf = 1e5, Δx = 1e4, Δy = 1e4)
        zonal_ocean = Subzero.Ocean(; grid, u = 1.0, v = 0.0, temp = 0.0)
        zero_atmos = Subzero.Atmos(; grid, u = 0.0, v = 0.0, temp = -20.0)
        domain = Subzero.Domain(;
            north = CollisionBoundary(North; grid),
            south = CollisionBoundary(South; grid),
            east = CollisionBoundary(East; grid),
            west = CollisionBoundary(West; grid),
        )
        floe = Subzero.Floe(
            [[
                [-1.75e4, 5e4],
                [-1.75e4, 7e4],
                [-1.25e4, 7e4], 
                [-1.25e4, 5e4],
                [-1.75e4, 5e4],
            ]],
            0.25,
            0.0,
        )
        area = floe.area
    # Standard Monte Carlo points for below floe - used to compare with MATLAB
        jldopen("inputs/test_mc_points.jld2", "r") do f
            floe.x_subfloe_points = f["X"]
            floe.y_subfloe_points = f["Y"]
        end
        
        modulus = 1.5e3*(sqrt(area) + sqrt(area))
        consts = Constants(E = modulus)
        floe_settings = FloeSettings()

    # stationary floe, uniform zonal ocean flow
        model1 = Model(
            grid,
            zonal_ocean,
            zero_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model1,
            10,
            consts,
            CouplingSettings(Δd = 2),
            floe_settings,
        )
        @test isapprox(model1.floes[1].fxOA/area, 2.9760, atol = 1e-3)
        @test isapprox(model1.floes[1].fyOA/area, 0.8296, atol = 1e-3)
        @test isapprox(model1.floes[1].trqOA/area, -523.9212, atol = 1e-3)

    # stationary floe, uniform meridional ocean flow
        meridional_ocean = Subzero.Ocean(; grid, u = 0.0, v = 1.0, temp = 0.0)
        model2 = Subzero.Model(
            grid,
            meridional_ocean,
            zero_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model2,
            10,
            consts,
            CouplingSettings(Δd = 2),
            floe_settings,
        )
        @test isapprox(model2.floes[1].fxOA/area, -0.8296, atol = 1e-3)
        @test isapprox(model2.floes[1].fyOA/area, 2.9760, atol = 1e-3)
        @test isapprox(model2.floes[1].trqOA/area, 239.3141, atol = 1e-3)

    # moving floe, uniform 0 ocean flow
        zero_ocean = Subzero.Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
        floe3 = deepcopy(floe)
        floe3.u = 0.25
        floe3.v = 0.1
        model3 = Subzero.Model(
            grid,
            zero_ocean,
            zero_atmos,
            domain,
            StructArray([floe3]),
        )
        Subzero.timestep_coupling!(
            model3,
            10,
            consts,
            CouplingSettings(Δd = 2),
            floe_settings,
        )
        @test isapprox(model3.floes[1].fxOA/area, -0.1756, atol = 1e-3)
        @test isapprox(model3.floes[1].fyOA/area, -0.1419, atol = 1e-3)
        @test isapprox(model3.floes[1].trqOA/area, 29.0465, atol = 1e-1)
        
        # stationary floe, diagonal atmos flow
        diagonal_atmos = Subzero.Atmos(; grid, u = -1, v = -0.5, temp = 0.0)
        model4 = Subzero.Model(
            grid,
            zero_ocean,
            diagonal_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model4,
            10,
            consts,
            CouplingSettings(Δd = 2),
            floe_settings,
        )
        @test isapprox(model4.floes[1].fxOA/area, -0.0013, atol = 1e-3)
        @test isapprox(model4.floes[1].fyOA/area, -6.7082e-4, atol = 1e-3)
        @test isapprox(model4.floes[1].trqOA/area, 0.2276, atol = 1e-3)

        # non-uniform ocean flow, zero atmos, stationary floe
        xgrid, ygrid = Subzero.grids_from_lines(
            grid.x0:grid.Δx:grid.xf,
            grid.y0:grid.Δy:grid.yf,
        )
        psi_ocn = 0.5e4*(sin.(4*(π/4e5).*xgrid) .* sin.(4*(π/4e5).*ygrid))
        non_unif_uocn = zeros(size(xgrid))
        non_unif_uocn[2:end, :] = -1e-4*(psi_ocn[2:end, :] .- psi_ocn[1:end-1, :])
        non_unif_vocn = zeros(size(ygrid))
        non_unif_vocn[:, 2:end] = 1e-4*(psi_ocn[:, 2:end] .- psi_ocn[:, 1:end-1])
        non_unif_ocean = Subzero.Ocean(;
            u = non_unif_uocn',
            v = non_unif_vocn',
            temp = zeros(size(xgrid)),
        )
        model5 = Subzero.Model(
            grid,
            non_unif_ocean,
            zero_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model5,
            10,
            consts,
            CouplingSettings(),
            floe_settings,
        )
        @test isapprox(model5.floes[1].fxOA/area, -0.0182, atol = 1e-3)
        @test isapprox(model5.floes[1].fyOA/area, 0.0392, atol = 1e-3)
        @test isapprox(model5.floes[1].trqOA/area, 23.6399, atol = 1e-3)


        # moving floe, non-uniform ocean, non-uniform atmos
        non_unif_atmos = Subzero.Atmos(;
            u = non_unif_uocn',
            v = non_unif_vocn',
            temp = zeros(size(xgrid)),
        )
        floe6 = deepcopy(floe)
        floe6.u = 0.5
        floe6.v = -0.5
        model6 = Subzero.Model(
            grid,
            non_unif_ocean,
            non_unif_atmos,
            domain,
            StructArray([floe6]),
        )
        Subzero.timestep_coupling!(
            model6,
            10,
            consts,
            CouplingSettings(),
            floe_settings,
        )
        @test isapprox(model6.floes[1].fxOA/area, -1.6300, atol = 1e-3)
        @test isapprox(model6.floes[1].fyOA/area, 1.1240, atol = 1e-3)
        @test isapprox(model6.floes[1].trqOA/area, 523.2361, atol = 2e-1)
    end
end