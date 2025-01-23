using Test, Subzero
import GeometryOps as GO
import GeometryOps.GeoInterface as GI

@testset "Directions" begin
    FT = Float64
    x0, xf, y0, yf = 0.0, 1e5, -5e4, 5e4
    Δx, Δy = (xf - x0) / 2, (yf - y0) / 2 
    # North boundary
    n_poly, n_val = Subzero._boundary_info_from_extent(North, FT, x0, xf, y0, yf)
    n_point_set = Set(GI.getpoint(n_poly))
    @test issetequal(n_point_set, Set(((x0 - Δx, yf), (x0 - Δx, yf + Δy), (xf + Δx, yf + Δy), (xf + Δx, yf))))
    @test n_val == yf
    # South boundary
    s_poly, s_val = Subzero._boundary_info_from_extent(South, FT, x0, xf, y0, yf)
    s_point_set = Set(GI.getpoint(s_poly))
    @test issetequal(s_point_set, Set(((x0 - Δx, y0), (x0 - Δx, y0 - Δy), (xf + Δx, y0 - Δy), (xf + Δx, y0))))
    @test s_val == y0
    # East boundary
    e_poly, e_val = Subzero._boundary_info_from_extent(East, FT, x0, xf, y0, yf)
    e_point_set = Set(GI.getpoint(e_poly))
    @test issetequal(e_point_set, Set(((xf, y0 - Δy), (xf, yf + Δy), (xf + Δx, y0 - Δy), (xf + Δx, yf + Δy))))
    @test e_val == xf
    # West boundary
    w_poly, w_val = Subzero._boundary_info_from_extent(West, FT, x0, xf, y0, yf)
    w_point_set = Set(GI.getpoint(w_poly))
    @test issetequal(w_point_set, Set(((x0, y0 - Δy), (x0, yf + Δy), (x0 - Δx, y0 - Δy), (x0 - Δx, yf + Δy))))
    @test w_val == x0
end

@testset "Boundaries" begin
    x0, xf, y0, yf = 0, 4e5, 0, 3e5
    grid = RegRectilinearGrid(; x0, xf, y0, yf, Δx = 1e4, Δy = 1e4)
    # Check Northern Boundary Values
    n_p = PeriodicBoundary(North; grid)
    n_o = OpenBoundary(North; x0, xf, y0, yf)
    # ensure val and poly fields are the same (and correct) regardless of constructor branch
    @test n_p.val == n_o.val == 3e5
    @test GO.equals(n_p.poly, n_o.poly)
    @test GO.area(n_p.poly) == 8e5 * 1.5e5
    north_points = Set(((-2e5, 3e5), (-2e5, 4.5e5), (6e5, 4.5e5), (6e5, 3e5)))
    @test issetequal(Set(GI.getpoint(n_p.poly)), north_points)
    # check type correctness
    @test typeof(n_p) == PeriodicBoundary{North, Float64}
    @test typeof(n_o) == OpenBoundary{North, Float64}

    # Check Eastern Boundary Values
    e_o = OpenBoundary(East; grid)
    e_c = CollisionBoundary(East; x0, xf, y0, yf)
    # ensure val and poly fields are the same (and correct) regardless of constructor branch
    @test e_o.val == e_c.val == 4e5
    @test GO.equals(e_o.poly, e_c.poly)
    @test GO.area(e_o.poly) == 2e5 * 6e5
    east_points = Set(((4e5, -1.5e5), (4e5, 4.5e5), (6e5, 4.5e5), (6e5, -1.5e5)))
    @test issetequal(Set(GI.getpoint(e_o.poly)), east_points)
    # check type correctness
    @test typeof(e_o) == OpenBoundary{East, Float64}
    @test typeof(e_c) == CollisionBoundary{East, Float64}
    
    # Check Southern Boundary Values
    s_p = PeriodicBoundary(South, Float64; grid)
    s_m = MovingBoundary(South, Float64; grid, u = 1.0, v = 2.0)
    # ensure val and poly fields are the same (and correct) regardless of constructor branch
    @test s_p.val == s_m.val == 0
    @test GO.equals(s_p.poly, s_m.poly)
    @test GO.area(s_p.poly) == 8e5 * 1.5e5
    south_points = Set(((-2e5, -1.5e5), (-2e5, 0.0), (6e5, 0.0), (6e5, -1.5e5)))
    s_p_points = Set(GI.getpoint(s_p.poly))
    @test issetequal(s_p_points, south_points)
    # check type correctness
    @test typeof(s_p) == PeriodicBoundary{South, Float64}
    @test typeof(s_m) == MovingBoundary{South, Float64}
    @test eltype(first(s_p_points)) == typeof(s_p.val) == Float64

    # Check Western Boundary Values
    w_c = CollisionBoundary(West, Float32; grid)
    w_m = MovingBoundary(West, Float32; x0, xf, y0, yf, u = 0.1, v = -0.1)
    # ensure val and poly fields are the same (and correct) regardless of constructor branch
    @test w_c.val == w_m.val == 0
    @test GO.equals(w_c.poly, w_m.poly)
    @test GO.area(w_c.poly, Float32) ≈ 2e5 * 6e5
    west_points = Set(((-2e5, -1.5e5), (-2e5, 4.5e5), (0.0, 4.5e5), (0.0, -1.5e5)))
    w_c_points = Set(GI.getpoint(w_c.poly))
    @test issetequal(w_c_points, west_points)
    # check type correctness
    @test typeof(w_c) == CollisionBoundary{West, Float32}
    @test typeof(w_m) == MovingBoundary{West, Float32}
    @test eltype(first(w_c_points)) == typeof(w_c.val) == Float32

    # Check _get_velocity boundary dispatch
    x, y = 0.0, 0.0  # (x, y) point does not matter as boundaries don't rotate
    @test Subzero._get_velocity(s_m, x, y) == (s_m.u, s_m.v) == (1.0, 2.0)
    @test Subzero._get_velocity(w_m, x, y) == (w_m.u, w_m.v) == (0.1f0, -0.1f0)
    @test Subzero._get_velocity(n_o, x, y) == (0.0, 0.0)
    @test Subzero._get_velocity(e_c, x, y) == (0.0, 0.0)
    @test Subzero._get_velocity(s_p, x, y) == (0.0, 0.0)

    # Check _update_boundary! boundary dispatch
    FT, Δt = Float64, 20
    # check that OpenBoundary does not move
    curr_val, curr_poly = n_o.val, n_o.poly
    Subzero._update_boundary!(n_o, Δt)
    @test n_o.val == curr_val && GO.equals(n_o.poly, curr_poly)
    # check that CollisionBoundary does not move
    curr_val, curr_poly = e_c.val, e_c.poly
    Subzero._update_boundary!(e_c, Δt)
    @test e_c.val == curr_val && GO.equals(e_c.poly, curr_poly)
    # check that PeriodicBoundary does not move
    curr_val, curr_poly = s_p.val, s_p.poly
    Subzero._update_boundary!(s_p, Δt)
    @test s_p.val == curr_val && GO.equals(s_p.poly, curr_poly)
    # check that Western MovingBoundary DOES move
    curr_val, curr_poly = w_m.val, w_m.poly
    Subzero._update_boundary!(w_m, Δt)
    Δx = Δt * w_m.u
    @test w_m.u != 0
    @test w_m.val == curr_val + Δx
    @test GO.equals(w_m.poly, Subzero._translate_poly(FT, curr_poly, Δx, 0.0))
    # check that Southern MovingBoundary DOES move
    curr_val, curr_poly = s_m.val, s_m.poly
    Subzero._update_boundary!(s_m, Δt)
    Δy = Δt * s_m.v
    @test s_m.v != 0
    @test s_m.val == curr_val + Δy
    @test GO.equals(s_m.poly, Subzero._translate_poly(FT, curr_poly, 0.0, Δy))

    # Check Periodic Compatability
    @test !Subzero._periodic_compat(n_p, s_m)
    @test !Subzero._periodic_compat(n_o, s_p)
    @test Subzero._periodic_compat(n_p, s_p)
    @test Subzero._periodic_compat(e_c, w_m)
end
