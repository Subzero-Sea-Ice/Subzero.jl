"""
    intersects(poly1::ArrayPolygon, poly2::ArrayPolygon)

Whether the polygons have a point in common, in their interiors or on their boundaries, as
`GO.intersects(poly1, poly2)`: the opposite of [`disjoint`](@ref GeometryOpsGPU.disjoint).
Polygons that only touch intersect. (For whether their interiors overlap, see
[`interiors_intersect`](@ref GeometryOpsGPU.interiors_intersect), for a rectangle.)
"""
intersects(poly1::ArrayPolygon, poly2::ArrayPolygon) = !disjoint(poly1, poly2)
