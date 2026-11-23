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

#=
Split the segment from `l_start` to `l_end` at every point where it meets `curve`, and report
whether the pieces run inside and/or outside the filled curve, as GeometryOps'
`_split_segment_interactions(Planar(), l_start, l_end, curve, in_curve, out_curve; exact)`:
`in_curve` and `out_curve` come in with the interactions found so far and go out with this
segment's added. It tests the midpoint of each piece.

GeometryOps collects the intersection points into a vector, sorts and deduplicates them, and
sorts them again by their distance from `l_start`. A kernel can't allocate that vector, so
the points are found one after another instead, in the same order
(`_next_segment_intersection`). That takes one pass over the curve per point, but the
segment usually meets the curve at only one to three points.
=#
function _split_segment_interactions(l_start::NTuple{2, T}, l_end::NTuple{2, T},
        curve::ArrayPolygon{<:AbstractArray{T}}, in_curve, out_curve; exact) where T
    p_start = l_start
    has_prev, key = false, (zero(T), l_start)
    while true
        found, next_key = _next_segment_intersection(l_start, l_end, curve, has_prev, key)
        p_end = found ? next_key[2] : l_end
        mid_val = _point_filled_curve_orientation((p_start .+ p_end) ./ 2, curve; exact)
        if mid_val == point_in
            in_curve = true
        elseif mid_val == point_out
            out_curve = true
        end
        found || break
        p_start, has_prev, key = p_end, true, next_key
    end
    return in_curve, out_curve
end

#=
The next point, after the one with key `prev` (or the first, if `has_prev` is false), where
the segment from `l_start` to `l_end` meets `curve`. The points are ordered by the key
`(distance from l_start, point)`, the order of GeometryOps' `_split_segment_interactions`:
it sorts the points (lexicographically) and removes duplicates, then sorts them by their
distance from `l_start`, keeping the order of points at the same distance. Returns whether
there is a next point, and its key. The points are those of
`GO.intersection_points(GI.Line([l_start, l_end]), curve)`, as in `intersection_points!`.
=#
function _next_segment_intersection(l_start::NTuple{2, T}, l_end::NTuple{2, T},
        curve::ArrayPolygon{<:AbstractArray{T}}, has_prev, prev) where T
    found, best = false, (zero(T), l_start)
    for j in 1:(npoint(curve) - 1)
        line_orient, (pt1, _), (pt2, _) = _intersection_point(T, (l_start, l_end),
            (getpoint(curve, j), getpoint(curve, j + 1)); exact = True())
        line_orient == line_out && continue
        found, best = _keep_next(found, best, (_euclid_distance(T, pt1, l_start), pt1),
            has_prev, prev)
        if line_orient == line_over
            found, best = _keep_next(found, best, (_euclid_distance(T, pt2, l_start), pt2),
                has_prev, prev)
        end
    end
    return found, best
end

# `key` if it comes after `prev` (when there is one) and before `best` (when there is one)
_keep_next(found, best, key, has_prev, prev) =
    (!has_prev || isless(prev, key)) && (!found || isless(key, best)) ? (true, key) :
        (found, best)

#=
Whether the line (a ring of an `ArrayPolygon`) interacts with the filled curve (the ring of
another) inside it, on its boundary, and outside it, as `(in_curve, on_curve, out_curve)`,
as `GO._line_filled_curve_interactions(Planar(), line, curve; exact, closed_line)`. With
`closed_line`, the line's last point connects to its first; a line whose last point repeats
its first is closed anyway.
=#
function _line_filled_curve_interactions(line::ArrayPolygon{<:AbstractArray{T}},
        curve::ArrayPolygon{<:AbstractArray{T}}; exact, closed_line = false) where T
    in_curve = false
    on_curve = false
    out_curve = false
    # the number of points of the line and the curve, without a repeated first point
    nl = npoint(line)
    nc = npoint(curve)
    first_last_equal_line = getpoint(line, 1) == getpoint(line, nl)
    first_last_equal_curve = getpoint(curve, 1) == getpoint(curve, nc)
    nl -= first_last_equal_line ? 1 : 0
    nc -= first_last_equal_curve ? 1 : 0
    closed_line |= first_last_equal_line
    # where the line's first point is
    l_start = getpoint(line, closed_line ? nl : 1)
    point_val = _point_filled_curve_orientation(l_start, curve; exact)
    if point_val == point_in
        in_curve = true
    elseif point_val == point_on
        on_curve = true
    else  # point_val == point_out
        out_curve = true
    end
    # where the line's segments meet the curve's
    for i in (closed_line ? 1 : 2):nl
        l_end = getpoint(line, i)
        c_start = getpoint(curve, nc)
        # every kind of interaction found already
        in_curve && on_curve && out_curve && break
        for j in 1:nc
            c_end = getpoint(curve, j)
            seg_val, _, _ = _seg_seg_orientation(l_start, l_end, c_start, c_end; exact)
            if seg_val != line_out
                # where the segments meet, the line is on the curve's boundary
                on_curve = true
                if seg_val == line_cross
                    # crossing the boundary, the line is both in and out of the curve
                    in_curve = true
                    out_curve = true
                else
                    if seg_val == line_over
                        sp = _point_segment_orientation(l_start, c_start, c_end)
                        lp = _point_segment_orientation(l_end, c_start, c_end)
                        if sp != point_in || lp != point_in
                            # the line goes over an end of the curve's segment, a hinge with
                            # the next segment
                            seg_val = line_hinge
                        end
                    end
                    if seg_val == line_hinge
                        #= A hinge may pass through other segments too, so find where the
                        pieces of the line's segment between the points where it meets the
                        curve are =#
                        in_curve, out_curve = _split_segment_interactions(
                            l_start, l_end, curve, in_curve, out_curve; exact,
                        )
                        # the whole segment is done
                        l_start = l_end
                        break
                    end
                end
            end
            c_start = c_end
        end
        l_start = l_end
    end
    return in_curve, on_curve, out_curve
end

#= Whether the line interacts with the polygon inside it, on its boundary, and outside it, as
`GO._line_polygon_interactions(Planar(), line, polygon; exact, closed_line)`. Polygons have
no holes, so that's the interactions with its ring. =#
_line_polygon_interactions(line::ArrayPolygon, polygon::ArrayPolygon; exact,
    closed_line = false) = _line_filled_curve_interactions(line, polygon; exact, closed_line)

#= Whether the two polygons interact as allowed and required, as GeometryOps'
`_polygon_polygon_process(Planar(), poly1, poly2; in_allow, on_allow, out_allow,
in_require, on_require, out_require, exact)`. If the polygons' bounding boxes don't meet,
that decides it. =#
@inline function _polygon_polygon_process(poly1::ArrayPolygon, poly2::ArrayPolygon; kw...)
    skip, returnval = _maybe_skip_disjoint_extents(poly1, poly2; kw...)
    skip && return returnval
    return _inner_polygon_polygon_process(poly1, poly2; kw...)
end

#= Whether the polygons' bounding boxes are disjoint, and if they are, the result, as
GeometryOps' `_maybe_skip_disjoint_extents(Planar(), a, b; ...)`. Returns
`(skip, returnval)`. =#
@inline function _maybe_skip_disjoint_extents(a::ArrayPolygon, b::ArrayPolygon;
        in_allow, on_allow, out_allow, in_require, on_require, out_require, kw...)
    ext_disjoint = Extents.disjoint(extent(a), extent(b))
    skip, returnval = if !ext_disjoint
        # can't tell anything about this case
        false, false
    elseif out_allow  # && ext_disjoint
        if in_require || on_require
            true, false
        else
            true, true
        end
    else  # !out_allow && ext_disjoint
        # points not allowed in exterior, but the polygons are disjoint
        true, false
    end
    return skip, returnval
end

#= `_polygon_polygon_process` once the bounding boxes meet, as GeometryOps'
`_inner_polygon_polygon_process`, without its steps for holes. =#
function _inner_polygon_polygon_process(poly1::ArrayPolygon, poly2::ArrayPolygon;
        in_allow, on_allow, out_allow, in_require, on_require, out_require, exact)
    in_req_met = !in_require
    on_req_met = !on_require
    out_req_met = !out_require
    # where the ring of poly1 is relative to poly2
    e1_in_p2, e1_on_p2, e1_out_p2 = _line_polygon_interactions(poly1, poly2;
        exact, closed_line = true)
    if e1_on_p2
        !on_allow && return false
        on_req_met = true
    end
    if e1_out_p2
        !out_allow && return false
        out_req_met = true
    end
    if !e1_in_p2
        # if the ring of poly1 isn't in poly2, whether poly1 surrounds poly2
        _, _, e2_out_e1 = _line_filled_curve_interactions(poly2, poly1;
            exact, closed_line = true)
        # if they really are disjoint, that's it
        e2_out_e1 && return in_req_met && on_req_met && out_req_met
    end
    # poly2 isn't outside of poly1, so the interiors of the polygons meet
    !in_allow && return false
    in_req_met = true
    return in_req_met && on_req_met && out_req_met
end
