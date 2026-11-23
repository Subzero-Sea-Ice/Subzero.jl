"""
    PointOrientation

Where a point is relative to a curve or polygon: `point_in`, `point_on` (on its boundary),
or `point_out`, as GeometryOps' `PointOrientation`.
"""
@enum PointOrientation point_in=1 point_on=2 point_out=3

#= Whether the point is in, on, or out of the polygon's ring, as
`GO._point_filled_curve_orientation(point, ring; exact)`: a ray-crossing test that also
detects points on an edge. `exact` is `False()` unless given, as for GeometryOps'
point-in-polygon tests; the clipping uses `True()`. =#
function _point_filled_curve_orientation(point::NTuple{2, T},
        p::ArrayPolygon{<:AbstractArray{T}}; exact = False()) where T
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
            f = Predicates.orient((xs, ys), (xe, ye), (x, y); exact)
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
function _point_polygon_process(point::NTuple{2, T}, p::ArrayPolygon{<:AbstractArray{T}};
        in_allow, on_allow, out_allow) where T
    x, y = point
    Extents.intersects(Extents.Extent(X = (x, x), Y = (y, y)), extent(p)) || return out_allow
    orientation = _point_filled_curve_orientation(point, p)
    orientation == point_out && return out_allow
    orientation == point_on && return on_allow
    return in_allow
end

#= Whether the point is `point_on` the segment from `start` to `stop` (one of its ends),
`point_in` it (any other of its points), or `point_out` of it, as
`GO._point_segment_orientation(Planar(), point, start, stop)`. =#
function _point_segment_orientation(point::NTuple{2, T}, start::NTuple{2, T},
        stop::NTuple{2, T}) where T
    x, y = point
    x1, y1 = start
    x2, y2 = stop
    Δx_seg = x2 - x1
    Δy_seg = y2 - y1
    Δx_pt = x - x1
    Δy_pt = y - y1
    if (Δx_pt == 0 && Δy_pt == 0) || (Δx_pt == Δx_seg && Δy_pt == Δy_seg)
        # the point is one of the segment's ends
        return point_on
    else
        # on the segment if parallel to it and between its ends
        on_line = _isparallel(Δx_seg, Δy_seg, Δx_pt, Δy_pt)
        !on_line && return point_out
        between_endpoints =
            (x2 > x1 ? x1 <= x <= x2 : x2 <= x <= x1) &&
            (y2 > y1 ? y1 <= y <= y2 : y2 <= y <= y1)
        !between_endpoints && return point_out
    end
    return point_in
end

#= The kind of intersection (a `LineOrientation`) of the segments from a1 to a2 and from b1
to b2, and the fractions α and β along them of its first point, as
`GO._seg_seg_orientation(Planar(), a1, a2, b1, b2; exact)`. Unlike GeometryOps, which
computes in Float64, it computes in the points' type. =#
function _seg_seg_orientation(a1::NTuple{2, T}, a2::NTuple{2, T}, b1::NTuple{2, T},
        b2::NTuple{2, T}; exact) where T
    seg_val, intr1, _ = _intersection_point(T, (a1, a2), (b1, b2); exact)
    (_, (α, β)) = intr1
    return seg_val, α, β
end
