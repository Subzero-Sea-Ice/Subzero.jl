using Test, Subzero
import GeometryOps as GO
import GeometryOps.GeoInterface as GI

@testset "Domain" begin
    grid = RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)
    north_periodic = PeriodicBoundary(North; grid)
    north_collision = CollisionBoundary(North; grid)
    south_periodic = PeriodicBoundary(South; grid)
    south_collision = CollisionBoundary(South; grid)
    east_periodic = PeriodicBoundary(East; grid)
    east_collision = CollisionBoundary(East; grid)
    west_periodic = PeriodicBoundary(West; grid)
    west_collision = CollisionBoundary(West; grid)
    
    # test basic periodic domain with no topography (creates default topography)
    domain1 = Domain(; north = north_periodic, south = south_periodic, east = east_periodic, west = west_periodic)
    @test typeof(domain1) <: Domain{Float64, <:AbstractBoundary, <:AbstractBoundary, <:AbstractBoundary, <:AbstractBoundary}
    @test domain1.north == north_periodic
    @test domain1.south == south_periodic
    @test domain1.east == east_periodic
    @test domain1.west == west_periodic
    @test typeof(domain1.topography) <: Subzero.TopographyField{Float64}
    @test isempty(domain1.topography)

    # test basic non-periodic domain with topography
    topography = initialize_topography_field(; coords = [[[(1e4, 1e4), (1e4, 2e4), (2e4, 2e4), (2e4, 1e4), (1e4, 1e4)]]])
    domain2 = Domain(; topography, north = north_collision, south = south_collision, east = east_collision, west = west_collision)
    @test domain2.north == north_collision
    @test domain2.south == south_collision
    @test domain2.east == east_collision
    @test domain2.west == west_collision
    @test typeof(domain2.topography) <: Subzero.TopographyField{Float64}
    @test length(domain2.topography) == 1
    @test GO.area(domain2.topography[1].poly) == 1e8

    # test error with wrong directions paired with wrong keywords --> error
    @test_throws MethodError Domain(; south = north_collision, north = south_collision, east = east_collision, west = west_collision)
    @test_throws MethodError Domain(; north = north_collision, south = south_collision, west = east_collision, east = west_collision)

    # domain with non-compatible periodic pairs  --> error
    @test_throws ArgumentError Domain(; north = north_collision, south = south_periodic, east = east_collision, west = west_collision)
    @test_throws ArgumentError Domain(; north = north_collision, south = south_collision, east = east_periodic, west = west_collision)

    # domain with north.val < south.val --> error
    north_low = PeriodicBoundary(North; x0 = 0.0, xf = 4e5, y0 = -4e5, yf = 0.0)
    south_high = PeriodicBoundary(South;  x0 = 0.0, xf = 4e5, y0 = 4e5, yf = 8e5)
    @test_throws ArgumentError Domain(; north = north_low, south = south_high, east = east_collision, west = west_collision)

    # domain with east < west --> error
    east_low = PeriodicBoundary(East; x0 = -4e5, xf = 0.0, y0 = 0.0, yf = 4e5)
    west_high = PeriodicBoundary(West;  x0 = 4e5, xf = 8e5, y0 = 0.0, yf = 4e5)
    @test_throws ArgumentError Domain(; north = north_collision, south = south_collision, east = east_low, west = west_high)

    # test Float32 Domain
    north_periodic_32 = PeriodicBoundary(North, Float32; grid)
    south_periodic_32 = PeriodicBoundary(South, Float32; grid)
    east_periodic_32 = PeriodicBoundary(East, Float32; grid)
    west_periodic_32 = PeriodicBoundary(West, Float32; grid)
    domain3 = Domain(; north = north_periodic_32, south = south_periodic_32, east = east_periodic_32, west = west_periodic_32)
    @test typeof(domain3) <: Domain{Float32, <:AbstractBoundary, <:AbstractBoundary, <:AbstractBoundary, <:AbstractBoundary}
    
    # domain with mixed types --> error
    @test_throws MethodError Domain(; north = north_periodic, south = south_periodic_32, east = east_periodic_32, west = west_periodic_32)
end