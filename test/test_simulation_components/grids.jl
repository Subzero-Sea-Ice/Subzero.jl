using Test, Subzero

@testset "CellFloes" begin
    # Empty CellFloes
    c1 = Subzero.CellFloes()
    @test isempty(c1.floeidx) && isempty(c1.Δx) && isempty(c1.Δy)
    @test typeof(c1) == Subzero.CellFloes{Float64}
    # Non-empty CellFloes
    c2 = Subzero.CellFloes(Float32; floeidx = [1, 2], Δx = [0.0, 0.0], Δy = [4.0, 4.0])
    @test length(c2.floeidx) == length(c2.Δx) == length(c2.Δy) == 2
    @test typeof(c2) == Subzero.CellFloes{Float32}
    # CellFloes with incorrect dimensions
    @test_throws AssertionError Subzero.CellFloes(Float32; floeidx = [1], Δx = [0.0, 0.0], Δy = [4.0, 4.0])
    # Empty CellFloes
    empty!(c2)
    @test isempty(c2.floeidx) && isempty(c2.Δx) && isempty(c2.Δy)
end

@testset "RegRectilinearGrid" begin
    # Non-square grid using constructor with Δx and Δy
    g1 = Subzero.RegRectilinearGrid(Float64; x0 = -10, xf = 10, y0 = -8, yf = 8, Δx = 2, Δy = 4)
    @test g1.Nx == 10
    @test g1.Ny == 4
    @test g1.x0 == -10
    @test g1.xf == 10
    @test g1.y0 == -8
    @test g1.yf == 8
    @test size(g1.floe_locations) == (11, 5)
    @test typeof(g1) == Subzero.RegRectilinearGrid{Float64}

    # Non-square grid using constructor with Nx and Ny
    g2 = Subzero.RegRectilinearGrid(Float32; x0 = -10, xf = 10, y0 = -8, yf = 8, Nx = 10, Ny = 4)
    @test g2.x0 == -10
    @test g2.xf == 10
    @test g2.y0 == -8
    @test g2.yf == 8
    @test g2.Δx == 2
    @test g2.Δy == 4
    @test size(g2.floe_locations) == (11, 5)
    @test typeof(g2) == Subzero.RegRectilinearGrid{Float32}

    # Test grid clipping to fit requested Δx and Δy when they divide evenly into (xf - x0) and (yf - y0)
    g3 = Subzero.RegRectilinearGrid(; x0 = 0, xf = 20, y0 = 0, yf = 16, Δx = 3, Δy = 3)
    @test g3.xf == 18
    @test g3.yf == 15
    @test g3.Nx == 6
    @test g3.Ny == 5
    @test size(g3.floe_locations) == (7, 6)
    @test typeof(g3) == Subzero.RegRectilinearGrid{Float64}
end