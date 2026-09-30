#=
Small version of the shear flow example: periodic boundaries (so ghost floes), collisions,
and coupling with a shear ocean current. Returns a function that creates the simulation for
a given backend. The floe field is generated once, because the Voronoi tesselation is not
reproducible with a fixed rng.
=#
function _shear_flow_simulation_factory(; FT = Float64)
    grid = RegRectilinearGrid(; x0 = 0.0, xf = 2e4, y0 = 0.0, yf = 2e4, Δx = 2e3, Δy = 2e3)
    domain = Domain(;
        north = PeriodicBoundary(North; grid),
        south = PeriodicBoundary(South; grid),
        east = PeriodicBoundary(East; grid),
        west = PeriodicBoundary(West; grid),
    )
    half_nx = fld(grid.Nx, 2)
    u_vec = [range(0, 0.5, length = half_nx); 0.5; range(0.5, 0, length = half_nx)]
    uvels = repeat(transpose(u_vec), outer = (grid.Ny + 1, 1))
    ocean = Ocean(; u = uvels, v = 0, temp = 0, grid)
    atmos = Atmos(; u = 0.0, v = 0.0, temp = -1.0, grid)
    floe_settings = FloeSettings(subfloe_point_generator = SubGridPointsGenerator(; grid, npoint_per_cell = 2))
    generator = VoronoiTesselationFieldGenerator(; nfloes = 10, concentrations = [0.75], hmean = 0.25, Δh = 0.0)
    floes = initialize_floe_field(FT; generator, domain, rng = Xoshiro(1), floe_settings)
    consts = Constants(; E = 1.5e3*(mean(sqrt.(floes.area)) + minimum(sqrt.(floes.area))))
    return function (backend)
        model = Model(; grid, ocean = deepcopy(ocean), atmos, domain, floes = deepcopy(floes))
        return Simulation(; model, consts, Δt = 20, nΔt = 100, floe_settings, rng = Xoshiro(1), backend)
    end
end

@testset "Simulation" begin
    make_simulation = _shear_flow_simulation_factory()
    @test make_simulation(Subzero.KernelAbstractions.CPU()).backend isa
        Subzero.KernelAbstractions.CPU
    gpu_backends = filter(b -> b isa Subzero.KernelAbstractions.GPU, test_backends())
    @testset "CPU and $backend give the same result" for backend in gpu_backends
        sims = [make_simulation(Subzero.KernelAbstractions.CPU()), make_simulation(backend)]
        with_logger(NullLogger()) do
            for sim in sims, tstep in 0:sim.nΔt
                timestep_sim!(sim, tstep)
            end
        end
        cpu_floes, gpu_floes = (sim.model.floes for sim in sims)
        @test length(cpu_floes) == length(gpu_floes)
        if length(cpu_floes) == length(gpu_floes)
            @test cpu_floes.id == gpu_floes.id
            @test reduce(vcat, cpu_floes.centroid) ≈ reduce(vcat, gpu_floes.centroid) rtol = 1e-8
            @test cpu_floes.u ≈ gpu_floes.u rtol = 1e-8
            @test cpu_floes.v ≈ gpu_floes.v rtol = 1e-8
            @test cpu_floes.ξ ≈ gpu_floes.ξ rtol = 1e-8
        end
    end
    @testset "run! does not change the global logger" begin
        grid = RegRectilinearGrid(; x0 = 0.0, xf = 1e4, y0 = 0.0, yf = 1e4, Δx = 2e3, Δy = 2e3)
        domain = Domain(;
            north = CollisionBoundary(North; grid),
            south = CollisionBoundary(South; grid),
            east = CollisionBoundary(East; grid),
            west = CollisionBoundary(West; grid),
        )
        ocean = Ocean(; u = 0.1, v = 0.0, temp = 0.0, grid)
        atmos = Atmos(; u = 0.0, v = 0.0, temp = 0.0, grid)
        floe_settings = FloeSettings()
        floes = initialize_floe_field(
            Float64,
            [[[[2e3, 2e3], [2e3, 4e3], [4e3, 4e3], [4e3, 2e3], [2e3, 2e3]]]],
            domain, 0.5, 0.0;
            floe_settings,
        )
        model = Model(; grid, ocean, atmos, domain, floes)
        sim = Simulation(; model, Δt = 10, nΔt = 5, floe_settings, name = "logger_test")
        logger_before = global_logger()
        # Default SubzeroLogger writes to a log file, which is closed afterwards
        run!(sim)
        @test global_logger() === logger_before
        @test isfile(joinpath("log", "logger_test.log"))
        # A given logger is used, but not closed
        io = IOBuffer()
        run!(sim; logger = SimpleLogger(io))
        @test global_logger() === logger_before
        @test isopen(io)
    end
end
