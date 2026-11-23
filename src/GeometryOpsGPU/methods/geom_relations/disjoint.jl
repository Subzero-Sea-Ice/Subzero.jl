# What disjoint allows, requires and uses, as GeometryOps' constants of the same names
const DISJOINT_ALLOWS = (in_allow = false, on_allow = false, out_allow = true)
const DISJOINT_REQUIRES = (in_require = false, on_require = false, out_require = false)
const DISJOINT_EXACT = (exact = False(),)

"""
    disjoint(poly1::ArrayPolygon, poly2::ArrayPolygon)

Whether the polygons have no point in common: neither their interiors nor their boundaries
meet, as `GO.disjoint(poly1, poly2)`. Like GeometryOps, it decides with inexact orientation
tests (`exact = False()`), which with Float32 polygons have Float32's tolerance. The polygons
must have the same float type.
"""
disjoint(poly1::ArrayPolygon{<:AbstractArray{T}},
    poly2::ArrayPolygon{<:AbstractArray{T}}) where T = _polygon_polygon_process(poly1, poly2;
        DISJOINT_ALLOWS..., DISJOINT_REQUIRES..., DISJOINT_EXACT...)
