"""
    PointOrientation

Where a point is relative to a curve or polygon: `point_in`, `point_on` (on its boundary),
or `point_out`, as GeometryOps' `PointOrientation`.
"""
@enum PointOrientation point_in=1 point_on=2 point_out=3

#= Whether the point is in, on, or out of the polygon's ring, as
`GO._point_filled_curve_orientation(point, ring; exact = False())`: a ray-crossing test
that also detects points on an edge. =#
function _point_filled_curve_orientation(point, p::ArrayPolygon)
    x, y = point
    n = npoint(p)
    n -= getpoint(p, 1) == getpoint(p, n) ? 1 : 0
    k = 0  # number of ray crossings
    xs, ys = getpoint(p, n)
    for j in 1:n
        xe, ye = getpoint(p, j)
        v1 = ys - y
        v2 = ye - y
        if !((v1 < 0 && v2 < 0) || (v1 > 0 && v2 > 0))
            u1, u2 = xs - x, xe - x
            f = Predicates.orient((xs, ys), (xe, ye), (x, y); exact = False())
            if v2 > 0 && v1 ≤ 0
                f == 0 && return point_on
                f > 0 && (k += 1)
            elseif v1 > 0 && v2 ≤ 0
                f == 0 && return point_on
                f < 0 && (k += 1)
            elseif v2 == 0 && v1 < 0
                f == 0 && return point_on
            elseif v1 == 0 && v2 < 0
                f == 0 && return point_on
            elseif v1 == 0 && v2 == 0
                u2 ≤ 0 && u1 ≥ 0 && return point_on
                u1 ≤ 0 && u2 ≥ 0 && return point_on
            end
        end
        xs, ys = xe, ye
    end
    return iseven(k) ? point_out : point_in
end

#= Whether the point's place relative to the polygon is allowed, as GeometryOps'
`_point_polygon_process(point, polygon; in_allow, on_allow, out_allow, exact = False())`:
a point outside the polygon's bounding box is outside, without the ray-crossing test. =#
function _point_polygon_process(point, p::ArrayPolygon; in_allow, on_allow, out_allow)
    x, y = point
    Extents.intersects(Extents.Extent(X = (x, x), Y = (y, y)), extent(p)) || return out_allow
    orientation = _point_filled_curve_orientation(point, p)
    orientation == point_out && return out_allow
    orientation == point_on && return on_allow
    return in_allow
end
