"""
    interiors_intersect(poly::ArrayPolygon, ext::Extents.Extent)

Whether the polygon overlaps the rectangle `ext` with a positive area, that is, whether
their interiors intersect. A polygon that only touches the rectangle doesn't. That's
`GO.intersects(poly, rect) && !GO.touches(poly, rect)`, and the same as
`!isempty(GO.intersection(poly, rect; target = GI.PolygonTrait()))`. GeometryOps has no
function of its own for it.

`ext` needs `X` and `Y` bounds, e.g. `GI.extent(rect)`, of the polygon's float type.

The interiors intersect if an edge of the polygon passes through the inside of the
rectangle. If no edge does, the inside of the rectangle is either entirely inside or
entirely outside the polygon, so they intersect if the polygon contains the rectangle's
centre.
"""
@inline function interiors_intersect(p::ArrayPolygon, ext::Extents.Extent)
    b = Extents.bounds(ext)
    return _interiors_intersect_rect(p, b.X[1], b.X[2], b.Y[1], b.Y[2])
end

# `interiors_intersect` for the rectangle xmin ≤ x ≤ xmax, ymin ≤ y ≤ ymax, whose bounds
# have the polygon's float type
@inline function _interiors_intersect_rect(p::ArrayPolygon{<:AbstractArray{T}},
        xmin::T, xmax::T, ymin::T, ymax::T) where T
    xc = (xmin + xmax) / 2
    yc = (ymin + ymax) / 2
    centre_inside = false
    # the ring is closed, so edge j goes from point j to point j + 1
    for j in 1:(npoint(p) - 1)
        x1, y1 = getpoint(p, j)
        x2, y2 = getpoint(p, j + 1)
        _edge_through_rectangle(x1, y1, x2, y2, xmin, xmax, ymin, ymax) && return true
        # even-odd rule: count the edges that a ray from the centre in +x direction crosses
        if (y1 > yc) != (y2 > yc) && xc < x1 + (yc - y1) * (x2 - x1) / (y2 - y1)
            centre_inside = !centre_inside
        end
    end
    return centre_inside
end

#=
Whether the edge from (x1, y1) to (x2, y2) passes through the inside of the rectangle
`xmin ≤ x ≤ xmax`, `ymin ≤ y ≤ ymax`, rather than only touching its sides or missing it.
Liang–Barsky clipping gives the part of the edge in the closed rectangle, as the
parameters t0 ≤ t ≤ t1 of the points (x1, y1) + t * (x2 - x1, y2 - y1). If the edge passes
through the inside, every point between the ends of that part is inside, so the midpoint
is. If it only touches the sides, the midpoint is on a side too.
=#
function _edge_through_rectangle(x1, y1, x2, y2, xmin, xmax, ymin, ymax)
    dx = x2 - x1
    dy = y2 - y1
    t0, t1 = zero(dx), one(dx)
    # the edge is on the inner side of each side of the rectangle where p * t ≤ q
    for (p, q) in ((-dx, x1 - xmin), (dx, xmax - x1), (-dy, y1 - ymin), (dy, ymax - y1))
        if iszero(p)  # parallel to this side
            q < 0 && return false
        elseif p < 0
            t0 = max(t0, q / p)
        else
            t1 = min(t1, q / p)
        end
    end
    t0 < t1 || return false
    tm = (t0 + t1) / 2
    xm = x1 + tm * dx
    ym = y1 + tm * dy
    return xmin < xm < xmax && ymin < ym < ymax
end
