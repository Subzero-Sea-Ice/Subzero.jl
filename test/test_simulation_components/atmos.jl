using Test, Subzero

@testset "Atmos" begin
    grid = Subzero.RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)
    Nx, Ny = grid.Nx, grid.Ny
    # Test Atmos with all constant real-input fields
    u_const, v_const, temp_const = 0.0, 0.1, -1.0
    atmos1 = Atmos(; grid, u = u_const, v = v_const, temp = temp_const)
    @test typeof(atmos1) == Atmos{Float64}
    @test size(atmos1.u) == size(atmos1.v) == size(atmos1.temp) == (Nx + 1, Ny + 1)
    @test all(atmos1.u .== u_const)
    @test all(atmos1.v .== v_const)
    @test all(atmos1.temp .== temp_const)

    # Test Atmos with all matrix-input fields
    seeded_rng = Xoshiro(1)
    u_matrix = rand(seeded_rng, Nx + 1, Ny + 1)
    v_matrix = rand(seeded_rng, Nx + 1, Ny + 1)
    temp_matrix = -1 .* rand(seeded_rng, Nx + 1, Ny + 1)
    atmos2 = Atmos(Float64; u = u_matrix, v = v_matrix, temp = temp_matrix)
    @test typeof(atmos2) == Atmos{Float64}
    @test size(atmos2.u) == size(atmos2.v) == size(atmos2.temp) == (Nx + 1, Ny + 1)
    @test all(atmos2.u .== u_matrix)
    @test all(atmos2.v .== v_matrix)
    @test all(atmos2.temp .== temp_matrix)
    
    # Test Atmos with mix of constant and matrix-input fields
    atmos3 = Atmos(Float32 ; grid, u = u_const, v = v_matrix, temp = temp_const)
    @test typeof(atmos3) == Atmos{Float32}
    @test size(atmos3.u) == size(atmos3.v) == size(atmos3.temp) == (Nx + 1, Ny + 1)
    @test all(atmos3.u .≈ u_const)
    @test all(atmos3.v .≈ v_matrix)
    @test all(atmos3.temp .≈ temp_const)

    # Test Atmos with const input fields but no grid
    @test_throws ArgumentError Atmos(; u = 0.0, v = 0.0, temp = 0.0)

    # Test field size constraints
    @test_throws ArgumentError Atmos(; u = zeros(3, 4), v = zeros(4, 3), temp = zeros(4, 3))
end
