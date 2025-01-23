using Test, Subzero
import GeometryOps as GO
import GeometryOps.GeoInterface as GI
import LibGEOS as LG

@testset "Topography" begin
    # Test basic polygon
    poly1 = GI.Polygon([[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]])
    topo1_Float64 = Subzero.TopographyElement(; poly = poly1)
    # check types
    @test typeof(topo1_Float64) == TopographyElement{Float64}
    @test GI.x(topo1_Float64.centroid) isa Float64 && GI.y(topo1_Float64.centroid) isa Float64
    @test topo1_Float64.rmax isa Float64
    # check values
    @test GO.equals(poly1, topo1_Float64.poly)
    @test all(topo1_Float64.centroid .≈ (0.5, 0.5))
    @test topo1_Float64.rmax ≈ sqrt(0.5^2 + 0.5^2)

    # Test basic polygon with Float32
    topo1_Float32 = Subzero.TopographyElement(Float32; poly = poly1)
    # check types
    @test typeof(topo1_Float32) == TopographyElement{Float32}
    @test GI.x(topo1_Float32.centroid) isa Float32 && GI.y(topo1_Float32.centroid) isa Float32
    @test topo1_Float32.rmax isa Float32

    # Test non-closed ring polygon
    poly2 = GI.Polygon([[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]])
    topo2 = Subzero.TopographyElement(; poly = poly2)
    # check close polygon
    @test GI.npoint(topo2.poly) == 5
    topo2_ring = GI.getexterior(topo2.poly)
    @test (GI.getpoint(topo2_ring, 1) == GI.getpoint(topo2_ring, 5))

    # Test polygon with hole
    poly3 = GI.Polygon([
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)],
        [(0.2, 0.2), (0.3, 0.2), (0.3, 0.3), (0.2, 0.3), (0.2, 0.2)],
    ])
    topo3 = Subzero.TopographyElement(; poly = poly3)
    # check polygon has no holes
    @test GI.nhole(topo3.poly) == 0

    # Test TopographyElement _get_velocity
    x, y = 0.5, 0.5
    @test Subzero._get_velocity(topo1_Float64, x, y) == (0.0, 0.0)

    # Test TopographyElement _normal_direction_correct!
    forces, points = ones(2, 2), ones(2, 2)
    Subzero._normal_direction_correct!(forces, points, topo1_Float64)
    @test forces == ones(2, 2)
end

@testset "TopographyField" begin
    coords1 = [[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]]
    coords2 = [[[0.5, 0.0], [0.5, 1.0], [1.5, 1.0], [1.5, 0.0], [0.5, 0.0]]]
    poly1 = GI.Polygon(coords1)
    poly2 = LG.Polygon(coords2)

    # Topography field of 2 identical polygons is reduced to one element
    topo_field_coords_64 = initialize_topography_field(; coords = [coords1, coords1])
    topo_field_polys_64 = initialize_topography_field(; polys = [poly1, poly1])
    # test field types
    @test typeof(topo_field_coords_64) == typeof(topo_field_polys_64)
    @test typeof(topo_field_coords_64) <: Subzero.TopographyField{Float64}
    # test polygons are the same and equal to poly1
    @test length(topo_field_coords_64) == length(topo_field_polys_64) == 1
    @test GO.equals(topo_field_coords_64.poly[1], poly1)
    @test GO.equals(topo_field_polys_64.poly[1], poly1)
    # test centroid and rmax
    @test topo_field_coords_64.centroid[1] == topo_field_polys_64.centroid[1] == [0.5, 0.5]
    @test topo_field_coords_64.rmax[1] == topo_field_polys_64.rmax[1] == sqrt(0.5^2 + 0.5^2)

    # Topography field of 2 overlapping polygons is make into two non-intersecting elements
    topo_field_coords_32 = initialize_topography_field(Float32; coords = [coords1, coords2])
    topo_field_polys_32 = initialize_topography_field(Float32; polys = [poly1, poly2])
    # test field types
    @test typeof(topo_field_coords_32) == typeof(topo_field_polys_32)
    @test typeof(topo_field_polys_32) <: Subzero.TopographyField{Float32}
    # test two constructor methods produce the same two polygons
    @test length(topo_field_coords_32) == length(topo_field_polys_32) == 2
    @test GO.equals(topo_field_coords_32.poly[1], topo_field_polys_32.poly[1])
    @test GO.equals(topo_field_coords_32.poly[2], topo_field_polys_32.poly[2])
    # test that the two polygons do not overlap
    (p1_x0, p1_xf), (p1_y0, p1_yf) = GI.extent(topo_field_coords_32.poly[1])
    @test (p1_x0, p1_xf) == (0.0, 0.5) && (p1_y0, p1_yf) == (0.0, 1.0)
    (p2_x0, p2_xf), (p2_y0, p2_yf) = GI.extent(topo_field_coords_32.poly[2])
    @test (p2_x0, p2_xf) == (0.5, 1.5) && (p2_y0, p2_yf) == (0.0, 1.0)
    # test centroid and rmax
    @test topo_field_coords_32.centroid[1] == topo_field_coords_32.centroid[1] == [0.25, 0.5]
    @test topo_field_coords_32.centroid[2] == topo_field_coords_32.centroid[2] == [1.0, 0.5]
    @test topo_field_coords_32.rmax[1] ≈ topo_field_coords_32.rmax[1] ≈ sqrt(0.25^2 + 0.5^2)
    @test topo_field_coords_32.rmax[2] ≈ topo_field_coords_32.rmax[2] ≈ sqrt(0.5^2 + 0.5^2)
end
