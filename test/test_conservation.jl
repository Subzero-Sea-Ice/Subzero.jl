function conservation_simulation(
    grid,
    domain,
    floes,
    smoothing = false,
    plot = false,
)
    Δt = 10
    ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
    atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = 0.0)
    model = Model(grid, ocean, atmos, domain, floes)
    dir = "output/conservation"
    initwriter = InitialStateOutputWriter(
        dir = dir,
        overwrite = true,
    )
    floewriter = FloeOutputWriter(
        10,
        dir = dir,
        overwrite = true,
    )
    writers = OutputWriters(
        initwriter,
        floewriter,
    )
    modulus = 1.5e3*(mean(sqrt.(floes.area)) + minimum(sqrt.(floes.area)))
    # No friction
    consts = Constants(E = modulus, μ = 0.0)
    # No ocean/atmosphere interactions
    coupling_settings = CouplingSettings(coupling_on = false)
    simplification_settings = SimplificationSettings(
        smooth_vertices_on = smoothing,
    )

    simulation = Simulation(
        model = model,
        consts = consts,
        Δt = 1,
        nΔt = 5000,
        verbose = false,
        coupling_settings = coupling_settings,
        simp_settings = simplification_settings,
        writers = writers,
    )
    run!(simulation)
    em_lists = check_energy_momentum_conservation_julia(
        joinpath(dir, "floes.jld2"),
        dir,
        plot
    )
    # Note that if any of these start with 0, this won't work
    return [vals[1] != 0 ?
        (vals[end] - vals[1])/vals[1] * 100 :
        NaN
        for vals in em_lists]
end


@testset "Conservation of Energy and Momentum" begin
    Δt = 10
    FT = Float64
    grid = RegRectilinearGrid(; x0 = -2e4, xf = 1e5, y0 = 0.0, yf = 1e5, Δx = 1e4, Δy = 1e4)
    collision_domain = Domain(; 
        north = CollisionBoundary(North; grid),
        south = CollisionBoundary(South; grid),
        east = CollisionBoundary(East; grid),
        west = CollisionBoundary(West; grid),
    )
    open_domain = Domain(;
        north = OpenBoundary(North; grid),
        south = OpenBoundary(South; grid),
        east = OpenBoundary(East; grid),
        west = OpenBoundary(West; grid),
    )
    topography = initialize_topography_field(; coords = [[[[-1e4, 0.0], [-2e4, 1e4], [-1e4, 1e4], [-1e4, 0.0]]]])
    open_domain_w_topography = Domain(;
        north = OpenBoundary(North; grid),
        south = OpenBoundary(South; grid),
        east = OpenBoundary(East; grid),
        west = OpenBoundary(West; grid),
        topography
    )
    rng = Xoshiro(1)
    floe1 = [[[2e4, 2e4], [2e4, 5e4], [5e4, 5e4], [5e4, 2e4], [2e4, 2e4]]]
    floe2 = [[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]
    floe3 = [[[5.5e4, 2e4], [5.25e4, 4e4], [5.75e4, 4e4], [5.5e4, 2e4]]]

    # Two blocks crashing head on - no rotation
    rng = Xoshiro(1)
    floe_settings = FloeSettings()
    head_on_floes = initialize_floe_field(
        FT,
        [floe1, floe2],
        open_domain, # Just affects shape, type doesn't matter
        0.25,
        0.0;
        rng = rng,
        floe_settings = floe_settings,
    )
    head_on_floes.u[1] = 0.15
    head_on_floes.u[2] = -0.1
    head_on_floes.v[1] = 0.02
    head_on_floes.v[2] = 0.02
    head_on_floes.ξ[1] = 1e-7
    # less than 1% change in both with open domain
    @test all(abs.(conservation_simulation(
        grid,
        open_domain,
        head_on_floes,
    )) .< 1)

    # Two blocks crashing offset - rotation
    rng = Xoshiro(1)
    offset_floes = initialize_floe_field(
        FT,
        [floe1, Subzero.translate(floe2, 0.0, 1e4)],
        open_domain, # Just affects shape, type doesn't matter
        0.25,
        0.0;
        rng = rng,
        floe_settings = floe_settings,
    )
    offset_floes.u[1] = 0.11
    offset_floes.u[2] = -0.1
    offset_floes.v[1] = 0.02
    offset_floes.v[2] = 0.02
    offset_floes.ξ[1] = 1e-7
    @test all(abs.(conservation_simulation(
        grid,
        open_domain,
        offset_floes,
    )) .< 1)

    # Two rectangular boxes with a triangle inbetween causing rotation
    rng = Xoshiro(1)
    rotating_floes = initialize_floe_field(
        FT,
        [floe1, floe2, floe3],
        open_domain, # Just affects shape, type doesn't matter
        0.25,
        0.0;
        rng = rng,
        floe_settings = floe_settings,
    )
    rotating_floes.u[1] = 0.11
    rotating_floes.u[2] = -0.1
    rotating_floes.ξ[3] = 1e-5
    rotating_floes.v[1] = 0.001
    rotating_floes.v[2] = 0.001
    rotating_floes.v[3] = 0.001
    @test all(abs.(conservation_simulation(
        grid,
        open_domain,
        rotating_floes,
    )) .< 1)

    # Three complex (many-sided, non-convex) floes hitting
    rng = Xoshiro(1)
    file = jldopen("inputs/floe_shapes.jld2", "r")
    complex_floes = initialize_floe_field(
        FT,
        [
            Subzero.translate(file["floe_vertices"][3], 0.0, 2e4),
            file["floe_vertices"][4],
            file["floe_vertices"][5],
        ],
        open_domain,
        0.25,
        0.0;
        rng = rng,
    )
    close(file)
    complex_floes.u[1] = 0.1
    complex_floes.v[2] = -0.2
    complex_floes.v[3] = 0.2
    # Slightly higher change in energy due to strage shapes
    @test all(abs.(
        conservation_simulation(
            grid,
            open_domain,
            complex_floes,)
    ) .< 2.1)

    # One non-convex block hits the wall and topography -> only check conservation of energy
    rng = Xoshiro(1)
    file = jldopen("inputs/floe_shapes.jld2", "r")
    floe_on_wall_topo = file["floe_vertices"][1]
    floe_on_wall_topo = Subzero.translate(floe_on_wall_topo, -1.75e4, -0.9e4)
    floe_arr = initialize_floe_field(
        FT,
        [floe_on_wall_topo],
        open_domain_w_topography,
        0.25,
        0.0;
        rng = rng,
    )
    close(file)
    floe_arr.u[1] = -0.09
    floe_arr.v[1] = -0.09
    @test abs(conservation_simulation(
        grid,
        open_domain_w_topography,
        floe_arr,
    )[1]) < 1
end