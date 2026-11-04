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
        origin_coords = translate_coords(
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
            MonteCarloPointsGenerator(Float32),
            GO.tuples(origin_poly, Float32),
            origin_centroid,
            area,
            Subzero.Status(),
            Xoshiro(1)
        )
        @test status3.tag == Subzero.active
        @test eltype(mc_x3) == eltype(mc_y3) == Float32

        # test generating sub-grid points for grid with Δx = Δy = 10
        point_generator = SubGridPointsGenerator(Float64; Δg = 10/sqrt(2))
        # Floe is smaller than grid cells --> centroid and vertices added
        square = Subzero.make_polygon([[
            [-2.5, -2.5],
            [-2.5, 2.5],
            [2.5, 2.5],
            [2.5, -2.5],
            [-2.5, -2.5],
        ]])
        square_centroid = GO.centroid(square)
        xpoints, ypoints, _ = Subzero.generate_subfloe_points(
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
        xpoints, ypoints, _ = Subzero.generate_subfloe_points(
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
        xpoints, ypoints, _ = Subzero.generate_subfloe_points(
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
        xpoints, ypoints, _ = Subzero.generate_subfloe_points(
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
         @test GO.equals(Subzero.center_cell_poly(Float64, 2, 3, grid, periodic_bound, periodic_bound),
            GI.Polygon([[[-9, -2], [-9, 2], [-7, 2], [-7, -2], [-9, -2]]]))
        @test GO.equals(Subzero.center_cell_poly(Float64, 1, 1, grid, open_bound, open_bound),
            GI.Polygon([[[-10, -8], [-10, -6], [-9, -6], [-9, -8], [-10, -8]]]))
        @test GO.equals(Subzero.center_cell_poly(Float64, 11, 6, grid, periodic_bound, periodic_bound),
            GI.Polygon([[[9, 10], [9, 14], [11, 14], [11, 10], [9, 10]]]))
        @test GO.equals(Subzero.center_cell_poly(Float64, 11, 6, grid, open_bound, open_bound),
            GI.Polygon([[[9, 8], [9, 8], [10, 8], [10, 8], [9, 8]]]))
        @test GO.equals(Subzero.center_cell_poly(Float64, 11, 6, grid, open_bound, periodic_bound),
            GI.Polygon([[[9, 8], [9, 8], [11, 8], [11, 8], [9, 8]]]))
        @test GO.equals(Subzero.center_cell_poly(Float64, 11, 6, grid, periodic_bound, open_bound),
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
        model1 = Model(;
            grid,
            ocean = zonal_ocean,
            atmos = zero_atmos,
            domain,
            floes = StructArray([deepcopy(floe)]),
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
        model2 = Subzero.Model(;
            grid,
            ocean = meridional_ocean,
            atmos = zero_atmos,
            domain,
            floes = StructArray([deepcopy(floe)]),
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
        model3 = Subzero.Model(;
            grid,
            ocean = zero_ocean,
            atmos = zero_atmos,
            domain,
            floes = StructArray([floe3]),
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
        model4 = Subzero.Model(;
            grid,
            ocean = zero_ocean,
            atmos = diagonal_atmos,
            domain,
            floes = StructArray([deepcopy(floe)]),
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
        model5 = Subzero.Model(;
            grid,
            ocean = non_unif_ocean,
            atmos = zero_atmos,
            domain,
            floes = StructArray([deepcopy(floe)]),
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
        model6 = Subzero.Model(;
            grid,
            ocean = non_unif_ocean,
            atmos = non_unif_atmos,
            domain,
            floes = StructArray([floe6]),
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

#=
CPU function that calculates subfloe point's cartesian coordiantes and index
within the grid and records them in pre-allocated input arrays, cart_vals 
and grid_idx. Returns the index of the last element in cart_vals and grid_idx
that corresponds to elements from given floe.
=#

function _ref_calc_subfloe_values!(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    grid,
    domain,
    cart_vals,
    grid_idx,
) where {FT<:AbstractFloat}
    # Translate/rotate monte carlo points to floe location/orientation
    α = floe.α
    j = 0  # index in output array
    for i in eachindex(floe.x_subfloe_points)  # index of monte carlo points
        px = cos(α)*floe.x_subfloe_points[i] -
            sin(α)*floe.y_subfloe_points[i]  # at origin
        py = sin(α)*floe.x_subfloe_points[i] +
            cos(α)*floe.y_subfloe_points[i]  # at origin
        x = px + floe.centroid[1]  # at centroid
        y = py + floe.centroid[2]  # at centroid
        # If point is in bounds, continue to find rest of values
        if in_bounds(x, y, grid, domain.north, domain.east)
            j += 1  # if added to outputs, move to next index in output array
            cart_vals[j, 1] = x
            cart_vals[j, 2] = y
            grid_idx[j, 1], grid_idx[j, 2] = find_center_cell_index(
                cart_vals[j, 1],
                cart_vals[j, 2],
                grid,
            )
        end
    end
    return j  # last spot in output arrays that has values for this floe
end

#=
    mc_interpolation(
        npoints,
        mc_grid_idx,
        grid,
        domain,
        atmos,
        ocean,
        coupling_settings,
    )

Create and returns interpolation objects for atmosphere u and v velocities, and
ocean u and v velocities, in addition to ocean's heatflux factor.
Inputs:
    npoints             <Int> number of monte carlo points to consider - the
                            number of rows to use in mc_cart and mc_grid_idx
    grid_idx            <Matrix{Int}> index of monte carlo points within the
                            grid - nx2 matrix of indices where the first column
                            is the grid column index and the second column is
                            the grid row index for cells centered on grid lines
    grid                <AbstractRectilinearGrid> model grid
    domain              <Domain> model domain
    atmos               <Atmos> model atmosphere
    ocean               <Ocean> model ocean
    coupling_settings   <CouplingSettings> simulation's coupling settings
Outputs:
    uatm_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the atompshere u velocity onto point
    vatm_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the atompshere v velocity onto point
    uocn_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean u velocity onto point
    vocn_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean v velocity onto point
    hflx_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean heatflux factor velocity onto point
=#
function mc_interpolation(
    npoints,
    grid_idx,
    grid,
    domain,
    atmos,
    ocean,
    coupling_settings,
)
    xidx = @view grid_idx[1:npoints, 1]
    yidx = @view grid_idx[1:npoints, 2]

    # Find knots and indices of knots for monte carlo interpolation
    xknots, xknot_idx = find_interp_knots(
        xidx,
        grid.Nx,
        grid.x0:grid.Δx:grid.xf,
        grid.xf - grid.x0,
        coupling_settings.Δd,
        domain.east,
    )
    yknots, yknot_idx = find_interp_knots(
        yidx,
        grid.Ny,
        grid.y0:grid.Δy:grid.yf,
        grid.yf - grid.y0,
        coupling_settings.Δd,
        domain.north,
    )

    knots = (xknots, yknots)

    # Atmos Interpolation objects for Monte Carlo Points
    uatm_interp = linear_interpolation(
        knots,
        @view(atmos.u[xknot_idx, yknot_idx]),
    )
    vatm_interp = linear_interpolation(
        knots,
        @view(atmos.v[xknot_idx, yknot_idx]),
    )
    
    # Ocean Interpolation objects for Monte Carlo Points
    uocn_interp = linear_interpolation(
        knots,
        @view(ocean.u[xknot_idx, yknot_idx]),
    )
    vocn_interp = linear_interpolation(
        knots,
        @view(ocean.v[xknot_idx, yknot_idx]),
    )
    hflx_interp = linear_interpolation(
        knots,
        @view(ocean.hflx_factor[xknot_idx, yknot_idx]),
    )

    return uatm_interp, vatm_interp, uocn_interp, vocn_interp, hflx_interp
end



#-------------- Ocean and Atmosphere on Ice --------------#
"""
    calc_atmosphere_forcing(...)

Calculates the stresses on a floe from the atmosphere above at given monte
carlo point.

## _Positional arguments_
- `xr::AbstractFloat`: x-coordiantes of points to interplate on
- `yr::AbstractFloat`: y-coordiantes of points to interplate on
- `upoint::AbstractFloat`: u velocity of floe at point
- `vpoint::AbstractFloat`: v velocity of floe at point
- `uatm_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the atompshere u velocity onto point
- `vatm_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the atompshere v velocity onto point
- $CONSTS_DEF

## _Returns_
- `τx_atm::AbstractFloat`: stress from atmosphere on floe in x-direction at given point
- `τy_atm::AbstractFloat`: stress from atmosphere on floe iny-direction at given point
"""
function calc_atmosphere_forcing(
    xr, 
    yr,
    upoint,
    vpoint,
    uatm_interp,
    vatm_interp,
    c,
)
    # Atmosphere velocities at monte carlo point
    uatm = uatm_interp(xr, yr) 
    vatm = vatm_interp(xr, yr) 

    # Stress on ice from atmopshere
    Δu_AI = uatm - upoint
    Δv_AI = vatm - vpoint
    norm = sqrt(Δu_AI^2 + Δv_AI^2)
    τx_atm = c.ρa * c.Cd_ia * norm * Δu_AI
    τy_atm = c.ρa * c.Cd_ia * norm * Δv_AI
    return τx_atm, τy_atm
end


"""
    calc_ocean_forcing!(...)

Calculates the stresses on a floe from the ocean above at given monte carlo
point.

## _Positional arguments_
- `xr::AbstractFloat`: x-coordiantes of points to interplate on
- `yr::AbstractFloat`: y-coordiantes of points to interplate on
- `upoint::AbstractFloat`: u velocity of floe at point
- `vpoint::AbstractFloat`: v velocity of floe at point
- `uocn_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the ocean u velocity onto point
- `vocn_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the ocean v velocity onto point
- `hflx_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the ocean heatflux factor velocity onto point
- `ma_ratio::AbstractFloat`: floe's mass to area ratio
- $CONSTS_DEF

## _Returns_
- `τx_ocn::AbstractFloat`: stress from ocean velocity on floe in x-direction at given point
- `τy_ocn::AbstractFloat`: stress from ocean velocity on floe in y-direction at given point
- `τx_pressure∇::AbstractFloat`: stress from ocean pressure gradient on floe in x-direction at given point
- `τy_pressure∇::AbstractFloat`: stress from ocean pressure gradient on floe in y-direction at given point
- `hflx_factor::AbstractFloat`: heatflux factor at given point from the heatflux factors of ocean below floe
"""
function calc_ocean_forcing!(
    xr,
    yr,
    upoint,
    vpoint,
    uocn_interp,
    vocn_interp,
    hflx_interp,
    ma_ratio,
    c,  # constants
)
    uocn = uocn_interp(xr, yr)
    vocn = vocn_interp(xr, yr)
    hflx_factor = hflx_interp(xr, yr)
    Δu_OI = uocn - upoint
    Δv_OI = vocn - vpoint
    norm = sqrt(Δu_OI^2 + Δv_OI^2)
    τx_ocn = c.ρo*c.Cd_io * norm * (cos(c.turnθ) * Δu_OI - sin(c.turnθ) * Δv_OI)
    τy_ocn = c.ρo*c.Cd_io * norm * (sin(c.turnθ) * Δu_OI + cos(c.turnθ) * Δv_OI)
    τx_pressure∇ = -ma_ratio * c.f * vocn
    τy_pressure∇ = ma_ratio * c.f * uocn
    return τx_ocn, τy_ocn, τx_pressure∇, τy_pressure∇, hflx_factor
end


"""
    floe_to_grid_info!(...)

Add force from the ice on ocean to ocean force fields (fx & fy) for each grid
cell and update ocean sea ice area fraction (si_area), representing total area
of sea ice in a given cell. Function is called for each sub-floe point.

## _Positional arguments_
- `floeidx::Int`: index of floe within model's floe array
- `xidx::Int`: grid x index that floe's point is within for grid centered on grid lines
- `yidx::Int`: grid column that floe's point is within for grid centered on grid lines
- `τx_ocn::AbstractFloat`: x-stress caused by ocean on point
- `τy_ocn::AbstractFloat`: y-stress caused by ocean on point
- $GRID_DEF
- $DOMAIN_DEF
- `cell_floes::Matrix{CellFloes}`: matrix of `CellFloes`, one for each grid cell
- `scells::Matrix{CellStresses}`: matrix of `CellStressess`, one for each grid cell
- `coupling_settings::CouplingSettings`: simulation's coupling settings

## _Returns_
- None. Updates `cell_floes` and `scells`.
"""
function floe_to_grid_info!(
    floeidx,
    xidx,
    yidx,
    τx_ocn::FT,
    τy_ocn::FT,
    grid,
    ns_bound,
    ew_bound,
    scells,
    coupling_settings,
) where {FT}
    # Determine grid cell point is in and if floe is shifted by periodic bounds
    shifted_xidx = shift_cell_idx(xidx, grid.Nx + 1, ew_bound)
    shifted_yidx = shift_cell_idx(yidx, grid.Ny + 1, ns_bound)
    Δx = (shifted_xidx - xidx) * (grid.Δx)
    Δy = (shifted_yidx - yidx) * (grid.Δy)
    if coupling_settings.two_way_coupling_on 
        # If two-way coupling, save stress on ocean per cell
        add_point!(
            grid.floe_locations[shifted_xidx, shifted_yidx],
            scells[shifted_xidx, shifted_yidx],
            floeidx,
            -τx_ocn,
            -τy_ocn,
            Δx,
            Δy,
        )
    else
        add_point!(
            grid.floe_locations[shifted_xidx, shifted_yidx],
            floeidx,
            Δx,
            Δy,
        )
    end
    return
end



"""
    calc_one_way_coupling!(...)

Preforms calculations needed for one way coupling by calculating floe's forcings
from ocean and atmosphere as well as the heatflux below a given floe.

## _Positional arguments_
- $FLOES_DEF
- $GRID_DEF
- `atmos::Atoms`: model's atmosphere
- `ocean::Ocean`: model's ocean
- $DOMAIN_DEF
- `coupling_settings::CouplingSettings`: simulation coupling settings
- $CONSTS_DEF

## _Returns_
- None. Update each floe's forces, torque, and heatflux factor from ocean/atmosphere.
"""
function calc_one_way_coupling!(
    floes::StructArray{Floe{FT}},
    grid,
    atmos,
    ocean,
    domain,
    coupling_settings,
    consts,
) where {FT}
    max_points = maximum(length, floes.x_subfloe_points)
    cart_vals = Matrix{FT}(undef, max_points, 2)
    grid_idx = Matrix{Int}(undef, max_points, 2)
    for i in eachindex(floes)
        # Monte carlo point cartesian coordinates and grid cell indices
        npoints = calc_subfloe_values!(
            get_floe(floes, i),
            grid,
            domain,
            cart_vals,
            grid_idx,
            validity_flag,
            npoints,
        )
        if npoints == 0
            floes.status[i].tag = remove
        else
            # Interpolaters for ocean and atmosphere
            uatm_int, vatm_int, uocn_int, vocn_int, hflx_int = mc_interpolation(
                npoints,
                grid_idx,
                grid,
                domain,
                atmos,
                ocean,
                coupling_settings,
            )

            # Add coriolis stress to total stress - same for every point
            xcoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.v[i]
            ycoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.u[i]
            tot_τx = npoints * xcoriolis
            tot_τy = -npoints * ycoriolis
            tot_τtrq = FT(0)
            tot_hflx_factor = FT(0)
            ma_ratio = floes.mass[i]/floes.area[i]
            # Determine total stress per-monte carlo point
            for j in 1:npoints
                # Monte carlo point properties
                xcentered = cart_vals[j, 1] - floes.centroid[i][1]
                ycentered = cart_vals[j, 2] - floes.centroid[i][2]
                θ = atan(ycentered, xcentered)
                rad = sqrt(xcentered^2 + ycentered^2)
                upoint = floes.u[i] - floes.ξ[i] * rad * sin(θ)
                vpoint = floes.v[i] + floes.ξ[i] * rad * cos(θ)
                # Stress at monte carlo point from ocean and atmosphere
                τx_atm, τy_atm = calc_atmosphere_forcing(
                    cart_vals[j, 1], 
                    cart_vals[j, 2],
                    upoint,
                    vpoint,
                    uatm_int,
                    vatm_int,
                    consts,
                )
                τx_ocn, τy_ocn, τx_p∇, τy_p∇, hflx_factor = calc_ocean_forcing!(
                    cart_vals[j, 1],
                    cart_vals[j, 2],
                    upoint,
                    vpoint,
                    uocn_int,
                    vocn_int,
                    hflx_int,
                    ma_ratio,
                    consts,
                )
                τx = τx_atm + τx_p∇ + τx_ocn
                τy = τy_atm + τy_p∇ + τy_ocn
                # Torque at monte carlo point
                τtrq = (-τx * sin(θ) + τy * cos(θ)) * rad
                # Add values to total stresses
                tot_τx += τx
                tot_τy += τy
                tot_τtrq += τtrq
                tot_hflx_factor += hflx_factor
                # Save floe info onto the grid
                floe_to_grid_info!(
                    i,
                    grid_idx[j, 1],
                    grid_idx[j, 2],
                    τx_ocn,
                    τy_ocn,
                    grid,
                    domain.north,
                    domain.east,
                    ocean.scells,
                    coupling_settings,
                )
            end
            # Average forces on ice floe
            floes.fxOA[i] = tot_τx/npoints * floes.area[i]
            floes.fyOA[i] = tot_τy/npoints * floes.area[i]
            floes.trqOA[i] = tot_τtrq/npoints * floes.area[i]
            floes.hflx_factor[i] = tot_hflx_factor/npoints
        end
    end
end
