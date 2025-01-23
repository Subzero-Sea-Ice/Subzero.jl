using Test, Subzero
import Random.Xoshiro

@testset "CellStresses" begin
    # Empty CellStresses
    c1 = Subzero.CellStresses()
    @test isempty(c1.τx) && isempty(c1.τy) && isempty(c1.npoints)
    @test typeof(c1) == Subzero.CellStresses{Float64}
    # Non-empty CellStresses
    c2 = Subzero.CellStresses(Float32; npoints = [1, 2], τx = [0.0, 0.0], τy = [4.0, 4.0])
    @test length(c2.npoints) == length(c2.τx) == length(c2.τy) == 2
    @test typeof(c2) == Subzero.CellStresses{Float32}
    # CellStresses with incorrect dimensions
    @test_throws AssertionError Subzero.CellStresses(Float32; npoints = [1], τx = [0.0, 0.0], τy = [4.0, 4.0])
    # Empty CellStresses
    empty!(c2)
    @test isempty(c2.npoints) && isempty(c2.τx) && isempty(c2.τy)
end

@testset "Ocean" begin
    grid = Subzero.RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)
    Nx, Ny = grid.Nx, grid.Ny
    # Test Ocean with all constant real-input fields
    u_const, v_const, temp_const = 0.0, 0.1, -1.0
    ocean1 = Ocean(; grid, u = u_const, v = v_const, temp = temp_const)
    @test typeof(ocean1) == Ocean{Float64}
    @test size(ocean1.u) == size(ocean1.v) == size(ocean1.temp) == (Nx + 1, Ny + 1)
    @test all(ocean1.u .== u_const)
    @test all(ocean1.v .== v_const)
    @test all(ocean1.temp .== temp_const)

    # make sure all other fields populated correctly
    zero_matrix = zeros(Nx + 1, Ny + 1)
    @test all(ocean1.hflx_factor .== zero_matrix)
    @test all(ocean1.τx .== zero_matrix)
    @test all(ocean1.τy .== zero_matrix)
    @test all(ocean1.si_frac .== zero_matrix)
    @test all(ocean1.dissolved .== zero_matrix)
    @test size(ocean1.scells) == (Nx + 1, Ny + 1)
    @test all([isempty(sc.τx) && isempty(sc.τy) && isempty(sc.npoints) for sc in ocean1.scells])

    # Test Ocean with all matrix-input fields
    seeded_rng = Xoshiro(1)
    u_matrix = rand(seeded_rng, Nx + 1, Ny + 1)
    v_matrix = rand(seeded_rng, Nx + 1, Ny + 1)
    temp_matrix = -1 .* rand(seeded_rng, Nx + 1, Ny + 1)
    ocean2 = Ocean(Float64; u = u_matrix, v = v_matrix, temp = temp_matrix)
    @test typeof(ocean2) == Ocean{Float64}
    @test size(ocean2.u) == size(ocean2.v) == size(ocean2.temp) == (Nx + 1, Ny + 1)
    @test all(ocean2.u .== u_matrix)
    @test all(ocean2.v .== v_matrix)
    @test all(ocean2.temp .== temp_matrix)
    
    # Test Ocean with mix of constant and matrix-input fields
    ocean3 = Ocean(Float32 ; grid, u = u_const, v = v_matrix, temp = temp_const)
    @test typeof(ocean3) == Ocean{Float32}
    @test size(ocean3.u) == size(ocean3.v) == size(ocean3.temp) == (Nx + 1, Ny + 1)
    @test all(ocean3.u .≈ u_const)
    @test all(ocean3.v .≈ v_matrix)
    @test all(ocean3.temp .≈ temp_const)

    # Test Ocean with const input fields but no grid
    @test_throws ArgumentError Ocean(; u = 0.0, v = 0.0, temp = 0.0)

    # Test ocean temperature constraint warning
    @test_logs (:warn, Subzero.OCEAN_TEMP_STR) Ocean(; grid, u = 0.0, v = 0.0, temp = 1.0)

    # Test field size constraints
    @test_throws ArgumentError Ocean(; u = zeros(3, 4), v = zeros(4, 3), temp = zeros(4, 3))
end
