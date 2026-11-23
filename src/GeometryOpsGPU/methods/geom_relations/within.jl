"""
    within(point, poly::ArrayPolygon)

Whether the point, a tuple `(x, y)`, is inside the polygon and not on its boundary, as
`GO.within(point, poly)`.
"""
within(point, p::ArrayPolygon) =
    _point_polygon_process(point, p; in_allow = true, on_allow = false, out_allow = false)
