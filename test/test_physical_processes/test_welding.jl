@testset "Bin floes" begin
    Δt = 10
    grid = RegRectilinearGrid(; x0 = 0.0, xf = 1e5, y0 = 0.0, yf = 1e5, Δx = 1e4, Δy = 1e4)

    open_domain = Subzero.Domain(;
        north = OpenBoundary(North; grid),
        south = OpenBoundary(South; grid),
        east = OpenBoundary(East; grid),
        west = OpenBoundary(West; grid),
    )
    periodic_domain = Subzero.Domain(;
        north = PeriodicBoundary(North; grid),
        south = PeriodicBoundary(South; grid),
        east = PeriodicBoundary(East; grid),
        west = PeriodicBoundary(West; grid),
    )
    half_open_periodic_domain = Subzero.Domain(;
        north = PeriodicBoundary(North; grid),
        south = PeriodicBoundary(South; grid),
        east = OpenBoundary(East; grid),
        west = OpenBoundary(West; grid),
    )

    coords = [
        [[[0.0, 1e4], [0.0, 4e4], [4e4, 4e4], [4e4, 1e4], [0.0, 1e4]]],   # Q4
        [[[1e4, 6e4], [1e4, 9e4], [4e4, 9e4], [4e4, 6e4], [1e4, 6e4]]],   # Q1
        [[[6e4, 6e4], [6e4, 9e4], [9e4, 9e4], [9e4, 6e4], [6e4, 6e4]]],   # Q2
        [[[6e4, 1e4], [6e4, 4e4], [9e4, 4e4], [9e4, 1e4], [6e4, 1e4]]],   # Q3
        [[[4e4, 4e4], [4e4, 6e4], [6e4, 6e4], [6e4, 4e4], [4e4, 4e4]]],   # Mid
        [[[9e4, 4e4], [9e4, 6e4], [11e4, 6e4], [11e4, 4e4], [9e4, 4e4]]], # Edge
        [[[4e4, -2e4], [4e4, 1e4], [6e4, 1e4], [6e4, -2e4], [4e4, -2e4]]] # Out
    ]
    logger = Logging.SimpleLogger(Logging.Error)
    floes = Logging.with_logger(logger) do  # suppress warning for floe outside domain
        initialize_floe_field(
            Float64,
            coords,
            periodic_domain,
            1.0,
            0.0;
        )
    end

    # Check assert with zero bins
    @test_throws AssertionError Subzero.bin_floe_centroids(
        floes,
        grid, open_domain,
        0, 1,  # Nx and Ny
    )
    @test_throws AssertionError Subzero.bin_floe_centroids(
        floes,
        grid, periodic_domain,
        1, 0,  # Nx and Ny
    )
    @test_throws AssertionError Subzero.bin_floe_centroids(
        floes,
        grid, half_open_periodic_domain,
        0, 0,  # Nx and Ny
    )

    # Check only one list with Nx = Ny = 1
    bins, nfloes = Subzero.bin_floe_centroids(
        floes,
        grid, periodic_domain,
        1, 1,  # Nx and Ny
    )
    @test size(bins) == (1, 1)
    @test length(bins[1, 1]) == 7
    @test nfloes[1, 1] == 7
    @test issetequal(Set(bins[1, 1]), Set(1:7)) 

    bins, nfloes = Subzero.bin_floe_centroids(
        floes,
        grid, open_domain,
        1, 1,  # Nx and Ny
    )
    @test size(bins) == (1, 1)
    @test length(bins[1, 1]) == 7
    @test nfloes[1, 1] == 6
    @test issetequal(Set(bins[1, 1]), Set(0:6)) 

    # Check floes sorted into two bins
    bins, nfloes = Subzero.bin_floe_centroids(
        floes,
        grid, open_domain,
        2, 1,  # Nx and Ny
    )
    @test size(bins) == (2, 1)
    @test length(bins[1, 1]) == 4
    @test nfloes[1, 1] == 2
    @test length(bins[2, 1]) == 4
    @test nfloes[2, 1] == 4
    @test issetequal(Set(bins[1, 1]), Set(0:2))
    @test issetequal(Set(bins[2, 1]), Set(3:6)) 

    bins, nfloes = Subzero.bin_floe_centroids(
        floes,
        grid, periodic_domain,
        2, 1,  # Nx and Ny
    )
    @test size(bins) == (2, 1)
    @test length(bins[1, 1]) == 4
    @test nfloes[1, 1] == 2
    @test length(bins[2, 1]) == 5
    @test nfloes[2, 1] == 5
    @test issetequal(Set(bins[1, 1]), Set(0:2))
    @test issetequal(Set(bins[2, 1]), Set(3:7)) 

    # Check floes sorted into 4 bins
    bins, nfloes = Subzero.bin_floe_centroids(
        floes,
        grid, half_open_periodic_domain,
        2, 2,  # Nx and Ny
    )
    @test size(bins) == (2, 2)
    @test length(bins[1, 1]) == 2
    @test nfloes[1, 1] == 1
    @test length(bins[1, 2]) == 2
    @test nfloes[1, 2] == 1
    @test length(bins[2, 1]) == 2
    @test nfloes[2, 1] == 2
    @test length(bins[2, 2]) == 3
    @test nfloes[2, 2] == 3
    @test issetequal(Set(bins[1, 1]), Set([0, 1]))
    @test issetequal(Set(bins[1, 2]), Set([0, 2]))
    @test issetequal(Set(bins[2, 1]), Set([4, 7]))
    @test issetequal(Set(bins[2, 2]), Set([3, 5, 6]))
end

@testset "Weld floes" begin
    Δt = 10
    grid = RegRectilinearGrid(; x0 = 0.0, xf = 1e5, y0 = 0.0, yf = 1e5, Δx = 1e4, Δy = 1e4)

    periodic_domain = Subzero.Domain(;
        north = OpenBoundary(North; grid),
        south = OpenBoundary(South; grid),
        east = OpenBoundary(East; grid),
        west = OpenBoundary(West; grid),
    )
    consts = Constants()
    floe_settings = FloeSettings()
    coupling_settings = CouplingSettings()
    coords = [
        [[[0.0, 0.0], [0.0, 5e4], [6e4, 5e4], [6e4, 0.0], [0.0, 0.0]]],
        [[[4e4, 0.0], [4e4, 5e4], [1e5, 5e4], [1e5, 0.0], [4e4, 0.0]]],
        [[[2e4, 4e4], [2e4, 8e4], [3e4, 8e4], [3e4, 4e4], [2e4, 4e4]]]
    ]
    floe_base = initialize_floe_field(
        Float64,
        coords,
        periodic_domain,
        1.0,
        0.0;
    )
    a1, a2, a3 = floe_base.area
    h1, h2, h3 = floe_base.height
    weld_settings = WeldSettings(
        weld_on = true,
        Δts = [700, 250, 100],
        Nxs = [1, 2, 1],
        Nys = [2, 2, 1],
        max_weld_area = 1e10,  # entire domain
        welding_coeff = 1000  # will weld for sure
    )
    # Test no floes welding because centroids in different bins (Nx = 2, Ny = 2)
    floes = deepcopy(floe_base)
    max_floe_id = Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        periodic_domain,
        weld_settings,
        floe_settings,
        2, # second set in weld setting lists (Nx = 2, Ny = 2)
        10,  # Δt
    )
    @test max_floe_id == 3
    @test Subzero.active == floes.status[1].tag ==
        floes.status[2].tag == floes.status[3].tag
    @test floes.area[1] == a1 && floes.area[2] == a2 && floes.area[3] == a3
    @test floes.height[1] == h1 && floes.height[2] == h2 &&
        floes.height[3] == h3

    # Test two floes welding with centroids in the same bin (Nx = 1, Ny = 2)
    floes = deepcopy(floe_base)
    max_floe_id = Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        periodic_domain,
        weld_settings,
        floe_settings,
        1, # first set in weld setting lists (Nx = 1, Ny = 2)
        10,  # Δt
    )
    @test max_floe_id == 4
    @test Subzero.active == floes.status[1].tag == floes.status[3].tag
    @test Subzero.remove == floes.status[2].tag
    @test floes.area[1] > a1
    @test floes.area[1] == 5e9
    @test floes.area[3] == a3
    @test floes.height[1] > h1
    @test floes.height[3] == h3
    @test floes.id[1] == 4
    @test length(floes.parent_ids[1]) == 2
    @test 1 in floes.parent_ids[1] && 2 in floes.parent_ids[1]

    # Test all three floes welding together (Nx = 1, Ny = 1)
    floes = deepcopy(floe_base)
    max_floe_id = Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        periodic_domain,
        weld_settings,
        floe_settings,
        3, # third set in weld setting lists (Nx = 1, Ny = 1)
        10,  # Δt
    )
    @test max_floe_id == 4
    @test Subzero.active == floes.status[1].tag
    @test Subzero.remove == floes.status[2].tag == floes.status[3].tag
    @test floes.area[1] == 5.3e9
    @test floes.area[1] > a1
    @test floes.height[1] > h1
    @test floes.id[1] == 4
    @test length(floes.parent_ids[1]) == 3
    @test 1 in floes.parent_ids[1] && 2 in floes.parent_ids[1] &&
        3 in floes.parent_ids[1]

    # All floes are too big to weld (use small maximum floe area for merging)
    floes = deepcopy(floe_base)
    max_floe_id = Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        periodic_domain,
        WeldSettings(
            weld_on = true,
            Δts = [100],
            Nxs = [1],
            Nys = [1],
            max_weld_area = 2.0e9,  # smaller than two largest floes
            welding_coeff = 1000  # will weld for sure
        ),
        floe_settings,
        1, # first set in weld setting lists (Nx = 1, Ny = 1)
        10,  # Δt
    )
    @test max_floe_id == 3
    @test Subzero.active == floes.status[1].tag ==
        floes.status[2].tag == floes.status[3].tag
    @test floes.area[1] == a1 && floes.area[2] == a2 && floes.area[3] == a3
    @test floes.height[1] == h1 && floes.height[2] == h2 &&
        floes.height[3] == h3

    # All floes are too small to weld
    floes = deepcopy(floe_base)
    max_floe_id = Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        periodic_domain,
        WeldSettings(
            weld_on = true,
            Δts = [100],
            Nxs = [1],
            Nys = [1],
            min_weld_area = 1e10,  # entire domain
            welding_coeff = 1000  # will weld for sure
        ),
        floe_settings,
        1, # first set in weld setting lists (Nx = 1, Ny = 1)
        10,  # Δt
    )
    @test max_floe_id == 3
    @test Subzero.active == floes.status[1].tag ==
        floes.status[2].tag == floes.status[3].tag
    @test floes.area[1] == a1 && floes.area[2] == a2 && floes.area[3] == a3
    @test floes.height[1] == h1 && floes.height[2] == h2 &&
        floes.height[3] == h3

    # Floe only welds with largest interaction (small max weld area)
    floes = deepcopy(floe_base)
    max_floe_id = Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        periodic_domain,
        WeldSettings(
            weld_on = true,
            Δts = [100],
            Nxs = [1],
            Nys = [1],
            max_weld_area = 5.1e9,  # smaller than all three floes together
            welding_coeff = 1000  # will weld for sure
        ),
        floe_settings,
        1, # first set in weld setting lists (Nx = 1, Ny = 1)
        10,  # Δt
    )
    @test max_floe_id == 4
    @test Subzero.active == floes.status[1].tag == floes.status[3].tag
    @test Subzero.remove == floes.status[2].tag
    @test floes.area[1] > a1
    @test floes.area[1] == 5e9
    @test floes.area[3] == a3
    @test floes.height[1] > h1
    @test floes.height[3] == h3
    @test length(floes.parent_ids[1]) == 2
    @test 1 in floes.parent_ids[1] && 2 in floes.parent_ids[1]
end