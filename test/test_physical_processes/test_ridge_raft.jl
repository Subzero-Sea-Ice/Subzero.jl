@testset "Ridging and Rafting" begin
    Δt = 10
    FT = Float64
    function update_height(floes, i, new_height, floe_settings)
        floes.height[i] = new_height
        floes.mass[i] = floes.area[i] * floes.height[i] * floe_settings.ρi
        floes.moment[i] = Subzero._calc_moment_inertia(
            Float64,
            Subzero.make_polygon(floes.coords[i]),
            floes.centroid[i],
            floes.height[i],
            ρi = floe_settings.ρi,
        )
        return
    end

    function assign_random_velocities!(floes)
        for i in eachindex(floes)
            floes.u[i] = (-1)^rand(0:1) * rand()
            floes.v[i] = (-1)^rand(0:1) * rand()
            floes.ξ[i] = (-1)^rand(0:1) * 0.05rand()
            floes.p_dxdt[i] = (-1)^rand(0:1) * rand()
            floes.p_dydt[i] = (-1)^rand(0:1) * rand()
            floes.p_dαdt[i] = (-1)^rand(0:1) * 0.05rand()
        end
        return
    end
    function setup_floes_with_inters(coords, domain, consts,
        collision_settings, lock, FT = Float64, Δx = nothing, Δy = nothing,
    )
        floes = initialize_floe_field(
            Float64,
            coords,
            domain,
            1.0,
            0.0;
        )
        if !isnothing(Δx)
            for i in eachindex(Δx)
                Subzero.translate!(floes.coords[i], Δx[i], Δy[i])
                floes.centroid[i][1] += Δx[i]
                floes.centroid[i][2] += Δy[i]
                floes.poly[i] = Subzero._translate_poly(FT, floes.poly[i], Δx[i], Δy[i])
            end
        end
        assign_random_velocities!(floes)
        add_ghosts!(floes, domain)
        Subzero.timestep_collisions!(  # Add interactions
            floes,
            length(floes),
            domain,
            consts,
            10,
            collision_settings, 
            lock,
        )
        # make sure all floes are active and not marked to fuse
        for i in eachindex(floes)
            floes.status[i].tag = Subzero.active
        end
        return floes
    end
    function calc_needed_momentum(floes)
        x_momentum, y_momentum = Subzero.calc_linear_momentum(
            floes.u,
            floes.v,
            floes.mass,
        )
        p_x_momentum, p_y_momentum = Subzero.calc_linear_momentum(
            floes.p_dxdt,
            floes.p_dydt,
            floes.mass,
        )
        return x_momentum, y_momentum, p_x_momentum, p_y_momentum
    end

    function conservation_of_momentum_tests(floes,
        x_momentum_init, y_momentum_init,
        p_x_momentum_init, p_y_momentum_init,
    )
        x_momentum_after, y_momentum_after,
        p_x_momentum_after, p_y_momentum_after = calc_needed_momentum(floes)
        
        @test abs((x_momentum_init - x_momentum_after) / x_momentum_init) < 1e-3
        @test abs((y_momentum_init - y_momentum_after) / y_momentum_init) < 1e-3
        @test abs((p_x_momentum_init - p_x_momentum_after) / p_x_momentum_init) < 1e-3
        @test abs((p_y_momentum_init - p_y_momentum_after) / p_y_momentum_init) < 1e-3
    end

    function test_floe_floe_rr_scenario(
        rr_settings,
        new_height1,
        new_height2,
        floe1_subsume,
        floe2_subsume,
        floes,
        domain,
        floe_settings,
        simp_settings,
    )
        @assert !floe1_subsume || !floe2_subsume "Both floe 1 and floe 2 can't subsume for testing"
        # Test floe 1 > min ridge height and floe 2 < min ridge height
        update_height(floes, 1, new_height1, floe_settings)
        update_height(floes, 2, new_height2, floe_settings)
        mass1, mass2 = floes.mass
        moment1, moment2 = floes.moment
        height1, height2 = floes.height
        cent1, cent2 = floes.centroid
        # initial momentum
        x_momentum_init, y_momentum_init,
        p_x_momentum_init, p_y_momentum_init = calc_needed_momentum(floes)
        # Ridge and raft floes
        Subzero.timestep_ridging_rafting!(
            floes,
            StructArray{Floe{Float64}}(undef, 0),
            domain,
            maximum(floes.id),
            rr_settings,
            floe_settings,
            simp_settings,
            10,
        )
        @test mass1 + mass2 == sum(floes.mass)
        conservation_of_momentum_tests(floes,
            x_momentum_init, y_momentum_init,
            p_x_momentum_init, p_y_momentum_init,
        )
        if floe1_subsume || floe2_subsume
            inter_polys =  Subzero.intersect_polys(Subzero.make_polygon(floes.coords[1]), Subzero.make_polygon(floes.coords[2]))
            @test sum(GO.area, inter_polys; init = 0.0) == 0  # floes DO NOT overlap anymore!
            @test floe1_subsume ?
                (mass1 < floes.mass[1] && mass2 > floes.mass[2]) :
                (mass1 > floes.mass[1] && mass2 < floes.mass[2])
            @test floe1_subsume ?
                (floes.height[1] > height1 && isapprox(floes.height[2], height2)) :
                (isapprox(floes.height[1], height1) && floes.height[2] > height2)
            @test floe1_subsume ?
                moment1 * floes.height[1] / height1 == floes.moment[1] :
                moment2 * floes.height[2] / height2 == floes.moment[2]
            @test floe1_subsume ?
                cent1 == floes.centroid[1] && cent2 != floes.centroid[2] :
                cent2 == floes.centroid[2] && cent1 != floes.centroid[1]
        else
            @test mass1 == floes.mass[1] && mass2 == floes.mass[2]
            @test floes.height[1] == height1 && floes.height[2] == height2
            @test moment1 == floes.moment[1] && moment2 == floes.moment[2]
            @test cent1 == floes.centroid[1] && cent2 == floes.centroid[2]
        end
    end
    function test_floe_domain_rr_scenario(rr_settings, does_raft, lose_mass, floes,
        domain, boundary_poly, topo_poly, bounds_overlap_area, topo_overlap_area,
        height1, height2, floe_settings, simp_settings)
        update_height(floes, 1, height1, floe_settings)
        update_height(floes, 2, height2, floe_settings)
        total_mass = sum(floes.mass)
        h1, h2 = floes.height
        area1, area2 = floes.area
        cent1, cent2 = floes.centroid

        x_momentum_init, y_momentum_init,
        p_x_momentum_init, p_y_momentum_init = calc_needed_momentum(floes)
        Subzero.timestep_ridging_rafting!(
            floes,
            StructArray{Floe{Float64}}(undef, 0),
            domain,
            maximum(floes.id),
            rr_settings,
            floe_settings,
            simp_settings,
            10,
        )
        conservation_of_momentum_tests(floes,
            x_momentum_init, y_momentum_init,
            p_x_momentum_init, p_y_momentum_init,
        )
        if does_raft
            if lose_mass
                @test total_mass > sum(floes.mass)
                @test isapprox(h1, floes.height[1]) && isapprox(h2, floes.height[2])
            else
                @test total_mass == sum(floes.mass)
                @test h1 < floes.height[1] && h2 < floes.height[2]
            end
            @test area1 - bounds_overlap_area == floes.area[1]
            @test area2 - topo_overlap_area == floes.area[2]
            @test cent1 != floes.centroid[1] && cent2 != floes.centroid[2]
            inter_polys = Subzero.intersect_polys(Subzero.make_polygon(floes.coords[1]), boundary_poly)
            @test sum(GO.area, inter_polys; init = 0.0) == 0
            inter_polys = Subzero.intersect_polys(Subzero.make_polygon(floes.coords[2]), topo_poly)
            @test sum(GO.area, inter_polys; init = 0.0) == 0
        else
            @test total_mass == sum(floes.mass)
            @test isapprox(h1, floes.height[1]) && isapprox(h2, floes.height[2])
            @test isapprox(area1,floes.area[1])
            @test isapprox(area2, floes.area[2])
            @test cent1 == floes.centroid[1] && cent2 == floes.centroid[2]
        end
    end

    function test_floe_ghost_rr_scenario(rr_settings, floes, height1, height2,
        floe1_subsume, domain, floe_settings, simp_settings,
    )
        update_height(floes, 1, height1, floe_settings)  # floe 1 will ridge onto floe 2
        update_height(floes, 3, height1, floe_settings)  # floe 1 has ghost floe 3
        update_height(floes, 2, height2, floe_settings)
        total_mass = floes.mass[1] + floes.mass[2]
        h1, h2, h3 = floes.height
        area1, area2, area3 = floes.area
        cent1, cent2, cent3 = floes.centroid
        x_momentum_init, y_momentum_init,
        p_x_momentum_init, p_y_momentum_init = calc_needed_momentum(floes[1:2])
        pieces_arr = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_arr,
            domain,
            maximum(floes.id),
            rr_settings,
            floe_settings,
            simp_settings,
            10,
        )
        # Test conservation of momentum
        conservation_of_momentum_tests(floes[1:2],
            x_momentum_init, y_momentum_init,
            p_x_momentum_init, p_y_momentum_init,
        )
        @test floes.u[1] == floes.u[3] && floes.v[1] == floes.v[3] &&
            floes.p_dxdt[1] == floes.p_dxdt[3] &&
            floes.p_dydt[1] == floes.p_dydt[3]

        @test length(floes) == 3
        # conservation of mass doesn't include ghost floes bc they're not real
        @test total_mass == floes.mass[1] + floes.mass[2]
        @test floe1_subsume ?
            h1 < floes.height[1] && h3 < floes.height[3] &&  isapprox(h1, h3) :
            isapprox(h1, floes.height[1]) && isapprox(h3, floes.height[3]) && isapprox(h1, h3)
        @test floe1_subsume ? 
            isapprox(h2, floes.height[2]) :
            h2 < floes.height[2]
        @test floe1_subsume ?
            isapprox(area1, floes.area[1]) && isapprox(area3, floes.area[3]) && isapprox(area1, area3) :
            area1 > floes.area[1] && area3 > floes.area[3] && isapprox(area1, area3)
        @test floe1_subsume ?
            area2 > floes.area[2] :
            isapprox(area2, floes.area[2])
        @test floe1_subsume ?
            cent2 != floes.centroid[2] :
            cent2 == floes.centroid[2]
        @test floe1_subsume ?
            cent1 == floes.centroid[1] && cent3 == floes.centroid[3] :
            cent1 != floes.centroid[1] && cent3 != floes.centroid[3]
        @test floes.centroid[1] .- floes.centroid[3] == cent1 .- cent3
    end

    # Setup for tests
    grid = RegRectilinearGrid(; x0 = 0, xf = 1e5, y0 = 0, yf = 1e5, Δx = 1e4, Δy = 1e4)
    topo_coords = [[[5e4, 5e4], [5e4, 7e4], [7e4, 7e4], [7e4, 5e4], [5e4, 5e4]]]
    collision_domain = Subzero.Domain(;
        north = CollisionBoundary(North; grid),
        south = CollisionBoundary(South; grid),
        east = CollisionBoundary(East; grid),
        west = CollisionBoundary(West; grid),
        topography = initialize_topography_field(; coords = [topo_coords])
    )
    periodic_domain = Subzero.Domain(;
        north = PeriodicBoundary(North; grid),
        south = PeriodicBoundary(South; grid),
        east = PeriodicBoundary(East; grid),
        west = PeriodicBoundary(West; grid),
    )
    boundary_poly = GO.UnionIntersectingPolygons()(Subzero.make_multipolygon([
        collision_domain.north.poly, collision_domain.south.poly,
        collision_domain.east.poly, collision_domain.west.poly]))
        
    topo_poly = Subzero.make_polygon(topo_coords)
    consts = Constants()
    floe_settings = FloeSettings(min_floe_area = 1e7)
    coupling_settings = CouplingSettings()
    simp_settings = SimplificationSettings()
    collision_settings = CollisionSettings(floe_floe_max_overlap = 1.0) # don't fuse
    lock = Threads.SpinLock()
    @testset "Floe-Floe Ridging and Rafting" begin
        coords = [
            [[[0.1e4, 0.1e4], [0.1e4, 2e4], [2e4, 2e4], [2e4, 0.1e4], [0.1e4, 0.1e4]]],
            [[[1.8e4, 1.8e4], [1.8e4, 4e4], [4e4, 4e4], [4e4, 1.8e4], [1.8e4, 1.8e4]]],
        ]
        floes_base = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock,
        )
        # Test scenario with no ridging or rafting
        no_rr_prob_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 0.0,  # no ridging
            raft_probability = 0.0,  # no rafting
        )
        floes = deepcopy(floes_base)
        # Test no ridging
        test_floe_floe_rr_scenario(
            no_rr_prob_settings,
            1.0,
            1.0,
            false,
            false,
            floes,
            collision_domain,
            floe_settings,
            simp_settings,
        )
        no_rr_frac_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,
            raft_probability = 1.0,
            min_overlap_frac = 1.0,  # need 100% overlap to ridge/raft
        )
        test_floe_floe_rr_scenario(
            no_rr_frac_settings,
            1.0,
            1.0,
            false,
            false,
            floes,
            collision_domain,
            floe_settings,
            simp_settings,
        )

        # # Ridging Tests
        ridge_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # force ridging
            raft_probability = 0.0,
            min_overlap_frac = 0.001,
        )
        floes = deepcopy(floes_base)
        # Test floe 2 ridging on top of floe 1
        test_floe_floe_rr_scenario(
            ridge_settings,
            1.0,
            0.1,
            true,  # floe 1 should subsume floe 2
            false,
            floes,
            collision_domain,
            floe_settings,
            simp_settings,
        )
        floes = deepcopy(floes_base)
        # Test floe 1 ridging on top of floe 2
        test_floe_floe_rr_scenario(
            ridge_settings,
            0.1,
            1.0,
            false,
            true,
            floes,
            collision_domain,
            floe_settings,
            simp_settings,
        )
        floes = deepcopy(floes_base)
        # Test both floes being too thin to ridge
        test_floe_floe_rr_scenario(
            ridge_settings,
            0.1,
            0.1,
            false,
            false,
            floes,
            collision_domain,
            floe_settings,
            simp_settings,
        )
        # Rafting Tests
        raft_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 0.0,  # force rafting
            raft_probability = 1.0,
            max_floe_raft_height = 1.0,
            min_overlap_frac = 0.001,
        )
        floes = deepcopy(floes_base)
        # Test floe 2 rafting on top of floe 1
        test_floe_floe_rr_scenario(
            raft_settings,
            1.0,
            0.001,  # bias so that floe 1 will subsume floe 2
            true,  # floe 1 should subsume floe 2
            false,
            floes,
            collision_domain,
            floe_settings,
            simp_settings,
        )
        # Test floe 2 rafting on top of floe 1
        floes = deepcopy(floes_base)
        test_floe_floe_rr_scenario(
            raft_settings,
            0.001, # bias so that floe 2 will subsume floe 1
            1.0,
            false,
            true,  # floe 2 should subsume floe 1
            floes,
            collision_domain,
            floe_settings,
            simp_settings,
        )
    end
    @testset "Floe-Domain Ridging and Rafting" begin
        coords = [
            [[[-0.1e4, -0.1e4], [-0.1e4, 2e4], [2e4, 2e4], [2e4, -0.1e4], [-0.1e4, -0.1e4]]],
            [[[3e4, 3e4], [3e4, 5e4], [5e4, 5e4], [5e4, 3e4], [3e4, 3e4]]]
        ]
        floes_base = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock, FT, [0.0, 0.5e4], [0.0, 0.5e4],
        )
        bounds_polys = Subzero.intersect_polys(Subzero.make_polygon(floes_base.coords[1]), boundary_poly)
        bounds_overlap_area = sum(GO.area, bounds_polys; init = 0.0)
        topos_polys = Subzero.intersect_polys(Subzero.make_polygon(floes_base.coords[2]), topo_poly)
        topo_overlap_area = sum(GO.area, topos_polys; init = 0.0)
 
        # Ridging with domain
        floes = deepcopy(floes_base)
        ridge_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # force ridging
            raft_probability = 0.0,
            min_overlap_frac = 0.001,
        )
        test_floe_domain_rr_scenario(ridge_settings, true, true, floes,
            collision_domain, boundary_poly, topo_poly, bounds_overlap_area,
            topo_overlap_area,  0.1, 0.1, floe_settings, simp_settings,
        )
        # Not ridging with domain
        floes = deepcopy(floes_base)
        test_floe_domain_rr_scenario(ridge_settings, false, true, floes,
            collision_domain, boundary_poly, topo_poly, bounds_overlap_area,
            topo_overlap_area, 2.0, 2.0, floe_settings, simp_settings,
        )
        # Rafting with domain
        floes = deepcopy(floes_base)
        rafting_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 0.0,  # force ridging
            raft_probability = 1.0,
            min_overlap_frac = 0.001,
        )
        test_floe_domain_rr_scenario(rafting_settings, true, true, floes,
            collision_domain, boundary_poly, topo_poly, bounds_overlap_area,
            topo_overlap_area, 0.1, 0.1, floe_settings, simp_settings,
        )
        # Not rafting with domain   
        floes = deepcopy(floes_base)
        test_floe_domain_rr_scenario(rafting_settings, false, true, floes,
            collision_domain, boundary_poly, topo_poly, bounds_overlap_area,
            topo_overlap_area, 0.3, 0.3, floe_settings, simp_settings,
        )

        # Ridging with domain where floe keeps mass
        ridge_keep_mass_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # force ridging
            raft_probability = 0.0,
            min_overlap_frac = 0.001,
            domain_gain_probability = 0.0,
        )
        floes = deepcopy(floes_base)
        test_floe_domain_rr_scenario(ridge_keep_mass_settings, true, false,
            floes, collision_domain, boundary_poly, topo_poly, bounds_overlap_area,
            topo_overlap_area, 1.0, 1.0, floe_settings, simp_settings,
        )

        coords = [[[[-0.1e4, 0.1e4], [-0.1e4, 9.9e4], [0.1e4, 9.9e4], [0.1e4, 0.1e4], [-0.1e4, 0.1e4]]]]
        floes_base = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock,
        )
        # Ridging with domain with small aspect ratio -> no piece saved
        floes = deepcopy(floes_base)
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            collision_domain,
            maximum(floes.id),
            ridge_settings,
            floe_settings,
            simp_settings,
            10,
        )
        @test floes.status[1].tag == Subzero.remove
        @test isempty(pieces_list)
    end

    @testset "Special Ridge Raft Cases" begin
        ridge_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # no ridging
            raft_probability = 0.0,  # no rafting
            min_overlap_frac = 0.0001,
        )
        coords = [
            [[[1e4, 1e4], [1e4, 5e4], [5e4, 5e4], [5e4, 1e4], [1e4, 1e4]]],
            [[[0.84e4, 2e4], [0.84e4, 2.5e4], [4e4, 2.5e4], [4e4, 2e4], [0.84e4, 2e4]]]
        ]
        floes_base = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock,
        )
        # floe 2 has a very small section outside of floe 1, too small to save
        floes = deepcopy(floes_base)
        update_height(floes, 2, 0.1, floe_settings)  # floe 2 will ridge onto floe 1
        total_mass = sum(floes.mass)
        h1, h2 = floes.height
        area1, area2 = floes.area
        cent1, cent2 = floes.centroid
        # Initial momentum 
        x_momentum_init, y_momentum_init,
        p_x_momentum_init, p_y_momentum_init = calc_needed_momentum(floes)
        # Run ridging and rafting
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            collision_domain,
            maximum(floes.id),
            ridge_settings,
            floe_settings,
            simp_settings,
            10,
        )
        conservation_of_momentum_tests(floes[1:1],
            x_momentum_init, y_momentum_init,
            p_x_momentum_init, p_y_momentum_init,
        )
        @test length(pieces_list) == 0
        @test isapprox(total_mass, floes.mass[1])
        @test floes.height[1] > h1
        @test floes.status[1].tag == Subzero.active
        @test floes.status[2].tag == Subzero.remove
        @test floes.centroid[1] == cent1

        # first floe overlaps with both other floes, and it will break when ridged with floe 2
        coords = [
            [[[2.75e4, 0.75e4], [0.75e4, 2.75e4], [1.25e4, 2.75e4], [3.25e4, 0.75e4], [2.75e4, 0.75e4]]],
            [[[0.1e4, 0.1e4], [0.1e4, 2.25e4], [2.25e4, 2.25e4], [2.25e4, 0.1e4], [0.1e4, 0.1e4]]],
            [[[2.5e4, 0.1e4], [2.5e4, 2.25e4], [3e4, 2.25e4], [3e4, 0.1e4], [2.5e4, 0.1e4]]],
        ]
        floes_base = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock,
        )
        # Test scenario with floe breaking into pieces
        floes = deepcopy(floes_base)
        update_height(floes, 1, 0.1, floe_settings)  # floe1 will ridge onto floe 2
        total_mass = sum(floes.mass)
        h1, h2, h3 = floes.height
        area1, area2, area3 = floes.area
        cent1, cent2, cent3 = floes.centroid
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            collision_domain,
            maximum(floes.id),
            ridge_settings,
            floe_settings,
            simp_settings,
            10,
        )
        @test length(pieces_list) == 1
        @test isapprox(total_mass, sum(floes.mass) + sum(pieces_list.mass))
        # Make sure floes 1 ridged onto floe 2 and broke
        @test isapprox(h1, floes.height[1]) && h2 < floes.height[2]
        @test isapprox(floes.height[1], pieces_list.height[1])
        @test cent2 == floes.centroid[2]
        # Make sure floe 3 wasn't ridged/rafted
        @test isapprox(h3, floes.height[3])
        @test cent3 == floes.centroid[3]
        # Make sure IDs are correct
        @test max_id == 5
        @test floes.id == [4, 2, 3] && pieces_list.id[1] == 5
        @test floes.parent_ids[1] == [1] && pieces_list.parent_ids[1] == [1]
        @test isempty(floes.parent_ids[2]) && isempty(floes.parent_ids[3])

        # Floe breaking on boundary and keeping mass
        ridge_keep_mass_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # force ridging
            raft_probability = 0.0,
            min_overlap_frac = 0.001,
            domain_gain_probability = 0.0,
        )
        coords = [
            [[[3e4, -0.2e4], [3e4, 0.2e4], [5e4, -0.1e4], [8e4, 0.2e4], [8e4, -0.2e4], [3e4, -0.2e4]]]
        ]
        base_floes = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock
        )
        floes = deepcopy(base_floes)
        total_mass = sum(floes.mass)
        h1 = floes.height[1]
        area1  = floes.area[1]
        cent1 = floes.centroid[1]
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            collision_domain,
            maximum(floes.id),
            ridge_keep_mass_settings,
            floe_settings,
            simp_settings,
            10,
        )
        @test length(pieces_list) == 1
        @test isapprox(total_mass, sum(floes.mass) + sum(pieces_list.mass))
        # Make sure floe ridged onto domain bounary and broke
        @test h1 < floes.height[1]
        @test h1 < pieces_list.height[1]
        @test isapprox(pieces_list.height[1], floes.height[1])
        @test cent1 != floes.centroid[1]
        # Make sure IDs are correct
        @test max_id == 3
        @test floes.id == [2] && pieces_list.id[1] == 3
        @test floes.parent_ids[1] == [1] && pieces_list.parent_ids[1] == [1]

        # Test parent ridging -> update ghost floes
        coords = [
            [[[-0.1e4, 0.1e4], [-0.1e4, 2e4], [2e4, 2e4], [2e4, 0.1e4], [-0.1e4, 0.1e4]]],
            [[[1.8e4, 1.8e4], [1.8e4, 4e4], [4e4, 4e4], [4e4, 1.8e4], [1.8e4, 1.8e4]]],
        ]
        base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
            collision_settings, lock
        )
        floes = deepcopy(base_floes)

        # Test parent-parent ridge, no breakage and update ghost
            # parent with ghost is subsumed by floe 2
        floes = deepcopy(base_floes)
        test_floe_ghost_rr_scenario(ridge_settings, floes, 0.1, 1.0,
            false, periodic_domain, floe_settings, simp_settings,
        )
            # parent with ghost subsumes floe 2 
        floes = deepcopy(base_floes)
        test_floe_ghost_rr_scenario(ridge_settings, floes, 1.0, 0.1,
            true, periodic_domain, floe_settings, simp_settings,
        )
        # Test parent-ghost ridge, no breakage and update parents
        coords = [
            [[[-0.1e4, 0.1e4], [-0.1e4, 2e4], [2e4, 2e4], [2e4, 0.1e4], [-0.1e4, 0.1e4]]],
            [[[8e4, 1.8e4], [8e4, 4e4], [9.92e4, 4e4], [9.92e4, 1.8e4], [8e4, 1.8e4]]],
        ]
        base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
            collision_settings, lock
        )
        # ghost (floe 3) is subsumed by floe 2
        floes = deepcopy(base_floes)
        test_floe_ghost_rr_scenario(ridge_settings, floes, 0.1, 1.0,
            false, periodic_domain, floe_settings, simp_settings,
        )
        # floe 2 is subsumed by ghost (floe 3)
        floes = deepcopy(base_floes)
        test_floe_ghost_rr_scenario(ridge_settings, floes, 1.0, 0.1,
            true, periodic_domain, floe_settings, simp_settings,
        )
    
        # Test parent-parent ridge, breakage
        coords = [
            [[[-1e4, 0.75e4], [0.75e4, 2.75e4], [1.25e4, 2.75e4], [-0.5e4, 0.75e4], [-1e4, 0.75e4]]],
            [[[0.1e4, 0.1e4], [0.1e4, 2.25e4], [2.25e4, 2.25e4], [2.25e4, 0.1e4], [0.1e4, 0.1e4]]],
            [[[9e4, 0.1e4], [9e4, 2.25e4], [9.5e4, 2.25e4], [9.5e4, 0.1e4], [9e4, 0.1e4]]],
        ]
        base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
            collision_settings, lock
        )
        floes = deepcopy(base_floes)
        update_height(floes, 1, 0.1, floe_settings)
        update_height(floes, 4, 0.1, floe_settings)  # ghost of floe 1
        total_mass = sum(floes.mass[1:3])
        h1, h2, h3, h4 = floes.height
        area1, area2, area3, area4 = floes.area
        cent1, cent2, cent3, cent4 = floes.centroid
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            periodic_domain,
            maximum(floes.id),
            ridge_settings,
            floe_settings,
            simp_settings,
            10,
        )
        # Make sure each piece is saved, ghosts are marked to remove and
        # don't interact with any other floes
        @test length(pieces_list) == 1
        @test isapprox(total_mass, sum(floes.mass[1:3]) + sum(pieces_list.mass))
        # Make sure floes 1 ridged onto floe 2 and broke
        @test isapprox(h1, floes.height[1], atol = 1e-15) && h2 < floes.height[2]
        @test isapprox(floes.height[1], pieces_list.height[1])
        @test isapprox(cent2, floes.centroid[2])
        # Make sure floe 3 wasn't ridged/rafted
        @test isapprox(h3, floes.height[3])
        @test isapprox(cent3, floes.centroid[3])
        # Make sure IDs are correct
        @test max_id == 5
        @test floes.id == [4, 2, 3, 1] && pieces_list.id[1] == 5
        @test floes.parent_ids[1] == [1] && pieces_list.parent_ids[1] == [1]
        @test isempty(floes.parent_ids[2]) && isempty(floes.parent_ids[3])
        # Make sure ghost floe isn't changed and is marked for remove_floe_overlap
        @test floes.status[4].tag == Subzero.remove
        @test isapprox(h4, floes.height[4])
        @test cent4 == floes.centroid[4]
        @test isapprox(area4, floes.area[4])

        # Test parent-ghost ridge, breakage
            # Make sure each piece is saved, ghosts are marked to remove and
            # don't interact with any other floes
        coords = [
            [[[8.5e4, 0.75e4], [10.75e4, 2.75e4], [11.25e4, 2.75e4], [9e4, 0.75e4], [8.5e4, 0.75e4]]],
            [[[8.5e4, 0.1e4], [8.5e4, 2.25e4], [9e4, 2.25e4], [9e4, 0.1e4], [8.5e4, 0.1e4]]],
            [[[0.1e4, 0.1e4], [0.1e4, 2.25e4], [2.25e4, 2.25e4], [2.25e4, 0.1e4], [0.1e4, 0.1e4]]],
        ]
        base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
            collision_settings, lock
        )
        floes = deepcopy(base_floes)
        update_height(floes, 1, 0.1, floe_settings)
        update_height(floes, 4, 0.1, floe_settings)  # ghost of floe 1
        total_mass = sum(floes.mass[1:3])
        h1, h2, h3, h4 = floes.height
        area1, area2, area3, area4 = floes.area
        cent1, cent2, cent3, cent4 = floes.centroid
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            periodic_domain,
            maximum(floes.id),
            ridge_settings,
            floe_settings,
            simp_settings,
            10,
        )
        # Make sure each piece is saved, ghosts are marked to remove and
        # don't interact with any other floes
        @test length(pieces_list) == 1
        @test isapprox(total_mass, sum(floes.mass[1:3]) + sum(pieces_list.mass))
        # Make sure floe 1 parent ridged onto floe 3 
        @test h3 < floes.height[3]
        @test isapprox(cent3, floes.centroid[3])
        # Make sure floes 1 ghost ridged onto floe 2 and broke
        @test isapprox(h1, floes.height[1], atol = 1e-15) && h2 < floes.height[2]
        @test isapprox(floes.height[1], pieces_list.height[1])
        @test cent2 == floes.centroid[2]
        # Make sure IDs are correct
        @test max_id == 5
        @test floes.id == [4, 2, 3, 1] && pieces_list.id[1] == 5
        @test floes.parent_ids[1] == [1] && pieces_list.parent_ids[1] == [1]
        @test isempty(floes.parent_ids[2]) && isempty(floes.parent_ids[3])
        # Make sure ghost floe is marked for remove_floe_overlap
        @test floes.status[4].tag == Subzero.remove
    end
end