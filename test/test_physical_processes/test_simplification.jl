@testset "Simplification" begin
    FT = Float64
    Δt = 10
    @testset "Dissolve Floes" begin
        grid = RegRectilinearGrid(; x0 = -1e5, xf = 1e5, y0 = 0.0, yf = 1e5, Δx = 1e4, Δy = 1e4)
        domain = Subzero.Domain(;
            north = CollisionBoundary(North; grid),
            south = CollisionBoundary(South; grid),
            east = PeriodicBoundary(East; grid),
            west = PeriodicBoundary(West; grid),
        )
        height = 0.25
        ρi = 920.0
        coords = [[
            [0.0, 5e4],
            [0.0, 8e4],
            [3e4, 8e4],
            [3e4, 5e4],
            [0.0, 5e4],
        ]]
        mass = 9e8 * height * ρi
        dissolved = zeros(Float64, 10, 20)  # 20x20 ocean grid
        # Add 2 floes in the middle of the grid -> masses added to cell
        floe = Floe(coords, height, 0.0,)
        Subzero.dissolve_floe!(floe, grid, domain, dissolved)
        @test dissolved[7, 12] == mass
        floe = Floe(Subzero.translate(coords, 2.5e3, 2.5e3), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolved)
        @test dissolved[7, 12] == 2mass
        # Add floe over periodic bound -> mass added to cell wrapped around grid
        floe = Floe(Subzero.translate(coords, 9e4, 0.0), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolved)
        @test dissolved[7, 1] == mass
        floe = Floe(Subzero.translate(coords, -1.2e5, 0.0), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolved)
        @test dissolved[7, 20] == mass
        total_mass = sum(dissolved)
        # Add floe over non-periodic bound -> mass not added since out of bounds
        floe = Floe(Subzero.translate(coords, 0.0, 6e4), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolved)
        @test total_mass == sum(dissolved)  # nothing was added
        floe = Floe(Subzero.translate(coords, 0.0, -7e4), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolved)
        @test total_mass == sum(dissolved)  # nothing was added
    end

    @testset "Fuse Floes" begin
        Δt = 10
        coords1 = [[
            [0.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [10.0, 0.0],
            [0.0, 0.0],
        ]]

        # Test two floes not intersecting -> will not fuse
        coords2 = deepcopy(coords1)
        Subzero.translate!(coords2, 20.0, 0.0)
        f1 = Floe(coords1, 0.5, 0.0)
        f2 = Floe(coords2, 0.5, 0.0)
        Subzero.fuse_two_floes!(
            f1,
            f2,
            10,
            FloeSettings(),
            2,
            Xoshiro(1),
        )
        @test f1.coords == coords1
        @test f2.coords == coords2

        # Test two floes intersecting -> will fuse into floe1 since same size
        Subzero.translate!(coords2, -13.0, 0.0)
        f2 = Floe(coords2, 0.75, 0.0)
        f1.id = 1
        f2.id = 2
        mass_tot = f1.mass + f2.mass
        f1.u = 0.1
        f1.v = 0.1
        f1.ξ = 0.1
        f1.p_dxdt = 0.002
        f1.p_dydt = 0.08
        f1.p_dαdt = 0.04
        f1.p_dudt = 0.01
        f2.p_dvdt = -0.007
        f2.p_dξdt = 0.01
        f2.u = 0.1
        f2.v = 0.2
        f2.ξ = 0
        f2.p_dudt = 0.02
        f2.p_dvdt = -0.005
        f2.p_dξdt = 0.05
        stress1_init = f1.stress_accum
        x_momentum_init, y_momentum_init = Subzero.calc_linear_momentum(
            [f1.u, f2.u],
            [f1.v, f2.v],
            [f1.mass, f2.mass],
        )
        spin_momentum_init, angular_momentum_init = Subzero.calc_angular_momentum(
            [f1.u, f2.u],
            [f1.v, f2.v],
            [f1.mass, f2.mass],
            [f1.ξ, f2.ξ],
            [f1.moment, f2.moment],
            [f1.centroid[1], f2.centroid[1]],
            [f1.centroid[2], f2.centroid[2]],
        )
        p_x_momentum_init, p_y_momentum_init = Subzero.calc_linear_momentum(
            [f1.p_dxdt, f2.p_dxdt],
            [f1.p_dydt, f2.p_dydt],
            [f1.mass, f2.mass],
        )
        p_spin_momentum_init, p_angular_momentum_init = Subzero.calc_angular_momentum(
            [f1.p_dxdt, f2.p_dxdt],
            [f1.p_dydt, f2.p_dydt],
            [f1.mass, f2.mass],
            [f1.p_dαdt, f2.p_dαdt],
            [f1.moment, f2.moment],
            [f1.centroid[1] - 10 * f1.p_dxdt, f2.centroid[1] - 10 * f2.p_dxdt],
            [f1.centroid[2] - 10 * f1.p_dydt, f2.centroid[2] - 10 * f2.p_dydt],
        )

        Subzero.fuse_two_floes!(
            f1,
            f2,
            10,
            FloeSettings(),
            4,
            Xoshiro(1),
        )
        @test f1.status.tag == Subzero.active
        @test f2.status.tag == Subzero.Subzero.remove
        @test f1.area == 170
        @test f1.mass == mass_tot  # conservation of mass
        @test f1.parent_ids == [1, 2]
        # conservation of momentum
        x_momentum_after, y_momentum_after = Subzero.calc_linear_momentum(
            [f1.u],
            [f1.v],
            [f1.mass],
        )
        spin_momentum_after, angular_momentum_after = Subzero.calc_angular_momentum(
            [f1.u],
            [f1.v],
            [f1.mass],
            [f1.ξ],
            [f1.moment],
            [f1.centroid[1]],
            [f1.centroid[2]],
        )
        p_x_momentum_after, p_y_momentum_after = Subzero.calc_linear_momentum(
            [f1.p_dxdt],
            [f1.p_dydt],
            [f1.mass],
        )
        p_spin_momentum_after, p_angular_momentum_after = Subzero.calc_angular_momentum(
            [f1.p_dxdt],
            [f1.p_dydt],
            [f1.mass],
            [f1.p_dαdt],
            [f1.moment],
            [f1.centroid[1] - 10 * f1.p_dxdt],
            [f1.centroid[2] - 10 * f1.p_dydt],
        )
        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-10)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-10)
        @test isapprox(p_x_momentum_init, p_x_momentum_after, atol = 1e-10)
        @test isapprox(p_y_momentum_init, p_y_momentum_after, atol = 1e-10)
        @test isapprox(
            spin_momentum_init + angular_momentum_init,
            spin_momentum_after + angular_momentum_after,
            atol = 1e-10,
        )
        @test isapprox(
            p_spin_momentum_init + p_angular_momentum_init,
            p_spin_momentum_after + p_angular_momentum_after,
            atol = 1e-10,
        )
        
        @test f1.stress_accum == (stress1_init * (f2.mass - mass_tot) .+ f2.stress_accum * f2.mass) / mass_tot

        # Test two floes intersecting -> will fuse into floe2 since bigger
        f1 = Floe(coords1, 0.5, 0.0)
        f3 = Floe(
            [[
                [0.0, 0.0],
                [0.0, 20.0],
                [20.0, 20.0],
                [20.0, 0.0],
                [0.0, 0.0],
            ]],
            0.55,
            0.0,
        )
        Subzero.fuse_two_floes!(
            f3,
            f1,
            10,
            FloeSettings(),
            3,
            Xoshiro(1),
        )
        @test f1.status.tag == Subzero.Subzero.remove
        @test f3.status.tag == Subzero.active

        # Test overall fuse floe functionality with set of 4 floes
        grid = RegRectilinearGrid(; x0 = -2.5e4, xf = 1e5, y0 = -2.5e4, yf = 1e5, Δx = 1e4, Δy = 1e4)
        open_domain_no_topo = Subzero.Domain(;
            north = OpenBoundary(North; grid),
            south = OpenBoundary(South; grid),
            east = OpenBoundary(East; grid),
            west = OpenBoundary(West; grid),
        )
        coords1 = [[  # large floe
            [0.0, 0.0],
            [0.0, 1e4],
            [1e4, 1e4],
            [1e4, 0.0],
            [0.0, 0.0],
        ]]
        coords2 = [[  # small floe
            [8e3, 5e3],
            [8e3, 8e3],
            [1.2e4, 8e3],
            [1.2e4, 5e3],
            [8e3, 5e3],
        ]]
        coords3 = [[  # large floe
            [1.1e4, 0.0],
            [1.1e4, 1e4],
            [2.1e4, 1e4],
            [2.1e4, 0.0],
            [1.1e4, 0.0]
        ]]
        coords4 = [[  # small floe
            [5e3, -2e3],
            [5e3, 3e3],
            [8e3, 3e3],
            [8e3, -2e3],
            [5e3, -2e3],
        ]]

        floe_arr = initialize_floe_field(
            FT,
            [coords1, coords2, coords3, coords4],
            open_domain_no_topo,
            0.5,
            0.0;
            floe_settings = FloeSettings(min_floe_area = 1e6),
            rng = Xoshiro(1),
        )
        floe_arr.status[1].tag = Subzero.fuse
        floe_arr.status[1].fuse_idx = [2, 4]
        floe_arr.status[2].tag = Subzero.fuse
        floe_arr.status[2].fuse_idx = [1, 3]
        floe_arr.status[3].tag = Subzero.fuse
        floe_arr.status[3].fuse_idx = [2]
        floe_arr.status[4].tag = Subzero.fuse
        floe_arr.status[4].fuse_idx = [4]

        max_floe_id = Subzero.fuse_floes!(
            floe_arr,
            4,
            FloeSettings(),
            10,
            Xoshiro(1),
        )
        floe3_area = floe_arr.area[3]
        @test max_floe_id == 6
        @test floe_arr.status[1].tag == Subzero.fuse
        @test floe_arr.id[1] == 6
        @test floe_arr.parent_ids[1] == [1, 2, 4]
        @test floe_arr.status[2].tag == Subzero.remove
        @test floe_arr.status[3].tag == Subzero.fuse
        @test floe_arr.status[4].tag == Subzero.remove
        @test floe_arr.area[3] == floe3_area  # small floes fused into floe 1
    end
    @testset "Smooth Floes" begin
        grid = RegRectilinearGrid(; x0 = -2.5e4, xf = 1e5, y0 = -2.5e4, yf = 1e5, Δx = 1e4, Δy = 1e4)
        open_domain_no_topo = Subzero.Domain(;
            north = OpenBoundary(North; grid),
            south = OpenBoundary(South; grid),
            east = OpenBoundary(East; grid),
            west = OpenBoundary(West; grid),
        )
        # Create complex floes
        file = jldopen("inputs/floe_shapes.jld2", "r")
        floe_coords = file["floe_vertices"][1:20]
        Subzero.translate!(floe_coords[2], 0.0, -1e3)
        floe_arr = initialize_floe_field(
            FT,
            floe_coords,
            open_domain_no_topo,
            0.5,
            0.0;
            floe_settings = FloeSettings(min_floe_area = 1e6),
            rng = Xoshiro(1),
        )
        close(file)

        floe_set1 = floe_arr[3:end]
        total_mass = sum(floe_set1.mass)
        nvertices = [length(c[1]) for c in floe_set1.coords]
        x_momentum_init, y_momentum_init = Subzero.calc_linear_momentum(
            floe_set1.u,
            floe_set1.v,
            floe_set1.mass,
        )
        spin_momentum_init, angular_momentum_init = Subzero.calc_angular_momentum(
            floe_set1.u,
            floe_set1.v,
            floe_set1.mass,
            floe_set1.ξ,
            floe_set1.moment,
            [c[1] for c in floe_set1.centroid],
            [c[2] for c in floe_set1.centroid],
        )
        p_x_momentum_init, p_y_momentum_init = Subzero.calc_linear_momentum(
            floe_set1.p_dxdt,
            floe_set1.p_dydt,
            floe_set1.mass,
        )
        p_spin_momentum_init, p_angular_momentum_init = Subzero.calc_angular_momentum(
            floe_set1.p_dxdt,
            floe_set1.p_dydt,
            floe_set1.mass,
            floe_set1.p_dαdt,
            floe_set1.moment,
            [f.centroid[1] - 10 * f.p_dxdt for f in floe_set1],
            [f.centroid[2] - 10 * f.p_dydt for f in floe_set1],
        )
        linear_e_init, rotation_e_init = Subzero.calc_kinetic_energy(
            floe_set1.u,
            floe_set1.v,
            floe_set1.mass,
            floe_set1.ξ,
            floe_set1.moment,
        )
        # smooth floes
        Subzero.smooth_floes!(
            floe_set1,
            open_domain_no_topo.topography,
            SimplificationSettings(max_vertices = 50),
            CollisionSettings(),
            FloeSettings(),
            10,
            Xoshiro(1),
        )

        for i in eachindex(floe_set1)
            # smooth floes if they have more than maximum number of vertices
            if nvertices[i] > 50
                @test length(floe_set1.coords[i][1]) < nvertices[i]
            else
                @test length(floe_set1.coords[i][1]) == nvertices[i]
            end
            @test floe_set1.status[i].tag == Subzero.active
        end
        # Confirm mass and momentum was conserved and energy wasn't gained
        @test total_mass == sum(floe_set1.mass)
        x_momentum_after, y_momentum_after = Subzero.calc_linear_momentum(
            floe_set1.u,
            floe_set1.v,
            floe_set1.mass,
        )
        spin_momentum_after, angular_momentum_after = Subzero.calc_angular_momentum(
            floe_set1.u,
            floe_set1.v,
            floe_set1.mass,
            floe_set1.ξ,
            floe_set1.moment,
            [c[1] for c in floe_set1.centroid],
            [c[2] for c in floe_set1.centroid],
        )
        p_x_momentum_after, p_y_momentum_after = Subzero.calc_linear_momentum(
            floe_set1.p_dxdt,
            floe_set1.p_dydt,
            floe_set1.mass,
        )
        p_spin_momentum_after, p_angular_momentum_after = Subzero.calc_angular_momentum(
            floe_set1.p_dxdt,
            floe_set1.p_dydt,
            floe_set1.mass,
            floe_set1.p_dαdt,
            floe_set1.moment,
            [f.centroid[1] - 10 * f.p_dxdt for f in floe_set1],
            [f.centroid[2] - 10 * f.p_dydt for f in floe_set1],
        )
        linear_e_after, rotation_e_after = Subzero.calc_kinetic_energy(
            floe_set1.u,
            floe_set1.v,
            floe_set1.mass,
            floe_set1.ξ,
            floe_set1.moment,
        )

        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-10)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-10)
        @test isapprox(p_x_momentum_init, p_x_momentum_after, atol = 1e-10)
        @test isapprox(p_y_momentum_init, p_y_momentum_after, atol = 1e-10)
        @test isapprox(
            spin_momentum_init + angular_momentum_init,
            spin_momentum_after + angular_momentum_after,
            atol = 1e-10,
        )
        @test isapprox(
            p_spin_momentum_init + p_angular_momentum_init,
            p_spin_momentum_after + p_angular_momentum_after,
            atol = 1e-10,
        )
        @test (linear_e_after + rotation_e_after) -
            (linear_e_init + rotation_e_init) <= 0

        # Two floes overlap, and one is cut into two pieces by topography
        floe_set2 = floe_arr[1:2]
        open_domain_with_topo = Subzero.Domain(;
            north = OpenBoundary(North; grid),
            south = OpenBoundary(South; grid),
            east = OpenBoundary(East; grid),
            west = OpenBoundary(West; grid),
            topography = initialize_topography_field(; coords = [[[
                [0.0, 1.05e4], [0.0, 1.15e4], [3e3, 1.15e4], [3e3, 1.05e4], [0.0, 1.05e4]]]])
        )
        og_f1_area = floe_set2.area[1]
        total_mass = sum(floe_set2.mass)
        nvertices = [length(c[1]) for c in floe_set2.coords]
        Subzero.smooth_floes!(
            floe_set2,
            open_domain_with_topo.topography,
            SimplificationSettings(max_vertices = 30),
            CollisionSettings(floe_floe_max_overlap = 0.05), # set low for test
            FloeSettings(),
            10,
            Xoshiro(1),
        )
        for i in eachindex(floe_set2)
            # smooth floes if they have more than maximum number of vertices
            @test length(floe_set2.coords[i][1]) < nvertices[i]
            @test floe_set2.status[i].tag == Subzero.fuse
        end
        # Test mass is conserved
        @test total_mass == sum(floe_set2.mass)
        # Test first floe was cut by topography and only larger piece was kept
        @test sum(GO.area, Subzero.intersect_polys(Subzero.make_polygon(floe_set2.coords[1]), open_domain_with_topo.topography.poly[1]); init = 0.0) == 0
        
        @test GO.area(Subzero.make_polygon(floe_set2.coords[1])) > 2og_f1_area/3
        # Test that both floes are tagged for fusion
        @test floe_set2.status[1].fuse_idx == [2]
        @test floe_set2.status[2].fuse_idx == [1]
    end
    @testset "Remove Floes" begin
        grid = RegRectilinearGrid(; x0 = -2.5e4, xf = 1e5, y0 = -2.5e4, yf = 1e5, Δx = 1e4, Δy = 1e4)

        open_domain_no_topo = Subzero.Domain(;
            north = OpenBoundary(North; grid),
            south = OpenBoundary(South; grid),
            east = OpenBoundary(East; grid),
            west = OpenBoundary(West; grid),
        )
        coords1 = [[  # large floe
            [0.0, 0.0],
            [0.0, 1e4],
            [1e4, 1e4],
            [1e4, 0.0],
            [0.0, 0.0],
        ]]
        coords2 = [[  # small floe
            [8e3, 5e3],
            [8e3, 8e3],
            [1.2e4, 8e3],
            [1.2e4, 5e3],
            [8e3, 5e3],
        ]]
        coords3 = [[  # large floe
            [1.1e4, 0.0],
            [1.1e4, 1e4],
            [2.1e4, 1e4],
            [2.1e4, 0.0],
            [1.1e4, 0.0]
        ]]
        coords4 = [[  # small floe
        [5e3, -2e3],
        [5e3, 3e3],
        [8e3, 3e3],
        [8e3, -2e3],
        [5e3, -2e3],
    ]]

        floe_arr = initialize_floe_field(
            FT,
            [coords1, coords2, coords3, coords4],
            open_domain_no_topo,
            0.5,
            0.0;
            floe_settings = FloeSettings(min_floe_area = 1e6),
            rng = Xoshiro(1),
        )
        floe_arr.height[4] = 0.05
        floe_arr.status[1].tag = Subzero.remove
        floe_arr.status[2].tag = Subzero.active
        floe_arr.status[3].tag = Subzero.fuse
        floe_arr.status[4].tag = Subzero.fuse
        floe_arr.id .= [1, 2, 3, 4]
        dissolve_mass = floe_arr.mass[2] + floe_arr.mass[4]

        dissolved = zeros(Float64, grid.Nx, grid.Ny)
        Subzero.remove_floes!(
            floe_arr,
            grid,
            open_domain_no_topo,
            dissolved,
            FloeSettings(
                min_floe_area = 1e8,
            )
        )
        @test length(floe_arr) == 1
        @test floe_arr.id == [3]
        @test floe_arr.status[1].tag == Subzero.active
        @test isempty(floe_arr.status[1].fuse_idx)
        @test sum(dissolved) == dissolve_mass
    end
end