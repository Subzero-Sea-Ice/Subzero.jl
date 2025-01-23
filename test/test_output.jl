
function test_basic_outputwriters()
    grid = RegRectilinearGrid(; x0 = -1e5, xf = 1e5, y0 = -1e5, yf = 1e5, Δx = 1e4, Δy = 1e4)
    ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
    atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = 0.0)
    domain = Domain(;
        north = OpenBoundary(North; grid),
        south = OpenBoundary(South; grid),
        east = OpenBoundary(East; grid),
        west = OpenBoundary(West; grid),
    )
    floe_coords = [[[7.5e4, 7.5e4], [7.5e4, 9.5e4], [9.5e4, 9.5e4], 
                    [9.5e4, 7.5e4], [7.5e4, 7.5e4]]]
    model = Model(
        grid,
        ocean,
        atmos,
        domain,
        StructArray([Floe(floe_coords, 0.5, 0.0)]),
    )

    dir = "output/sim"
    initwriter = InitialStateOutputWriter(
        dir = dir,
        filename = "initial_state.jld2",
        overwrite = true,
    )
    gridwriter = GridOutputWriter(
        100,
        grid,
        (10, 5),
        dir = dir,
        filename = "grid.nc",
        overwrite = true,
    )
    floewriter = FloeOutputWriter(
        50;
        outputs = [:status, :coords, :area, :mass, :u, :v],
        dir = dir,
        filename = "floe.jld2",
        overwrite = true,
    )
    checkpointwriter = CheckpointOutputWriter(
        250,
        dir = dir,
        overwrite = true,
    )
    writers = OutputWriters(
        initwriter,
        gridwriter,
        floewriter,
        checkpointwriter,
    )
    simulation = Simulation(
        model = model,
        consts = Constants(),
        nΔt = 500,
        writers = writers)

    run!(simulation)

    # Test Initial State Output Writer
    fn = joinpath(dir, "initial_state.jld2")
    file = jldopen(fn, "r")
    @test keys(file) == ["sim"]
    @test typeof(file["sim"]) <: Subzero.Simulation
    close(file)
    rm("output/sim/initial_state.jld2")

    # Test Checkpoint Output Writer
    fn = joinpath(dir, "checkpoint.jld2")
    file = jldopen(fn, "r")
    @test Set(keys(file)) == Set(["ocean", "metadata", "atmos", "floes"])
    @test typeof(file["ocean/0"]) <: Subzero.Ocean
    @test typeof(file["atmos/0"]) <: Subzero.Atmos
    @test eltype(file["floes/0"]) <: Subzero.Floe
    @test keys(file["ocean"]) == ["0", "250", "500"]
    close(file)
    rm(fn)

    # Test Floe Output Writer
    fn = joinpath(dir, "floe.jld2")
    file = jldopen(fn, "r")
    @test Set(keys(file)) ==
        Set(["status", "coords", "area", "mass", "u", "v", "metadata"])
    @test length(keys(file["status"])) == 11
    close(file)
    rm(fn)

    # Test Grid Output Writer
    fn = joinpath(dir, "grid.nc")
    file = Dataset(fn, "r")
    @test Set(keys(file)) ==
        Set(["u_grid", "v_grid", "dudt_grid", "dvdt_grid", "overarea_grid",
            "mass_grid", "area_grid", "height_grid", "si_frac_grid",
            "stress_xx_grid", "stress_yx_grid", "stress_xy_grid",
            "stress_yy_grid", "stress_eig_grid", "strain_ux_grid",
            "strain_vx_grid", "strain_uy_grid", "strain_vy_grid", "time", "x",
            "y"]
        )
    @test length(file["x"][:]) == 10
    @test length(file["y"][:]) == 5
    @test length(file["time"][:]) == 6
    close(file)
    rm(fn)
end

@testset "Outputwriters" begin
    test_basic_outputwriters()
end

