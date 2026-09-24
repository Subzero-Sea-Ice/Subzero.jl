@testset "Simulation" begin
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
