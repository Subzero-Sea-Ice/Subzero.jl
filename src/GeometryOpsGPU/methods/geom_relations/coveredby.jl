"""
    coveredby(point, poly::ArrayPolygon)

Whether the point, a tuple `(x, y)` of the polygon's float type, is inside the polygon or on
its boundary, as `GO.coveredby(point, poly)`.
"""
coveredby(point, p::ArrayPolygon) =
    _point_polygon_process(point, p; in_allow = true, on_allow = true, out_allow = false)
