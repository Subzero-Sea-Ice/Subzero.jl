#=
Geometry on polygons stored in plain arrays, so that it can run in GPU kernels: no
allocations and no GeometryOps or GeoInterface objects. Each function gives the same result
as the GeometryOps function named in its comment, bit for bit: it does the same operations
in the same order and in the same types, including where GeometryOps converts to Float64.
Rings are closed (the last point repeats the first) and have no holes, like floes.

Like GeometryOps, the orientation tests (`_orient_inexact`, `_orient_exact`) use Float64 even
for Float32 floes, so these functions can't run on GPUs without Float64 (Metal) yet. The
exact test is AdaptivePredicates' `orient2p`, which runs in kernels.
=#

#=
A closed ring of points stored in a 3D array: point `j` of ring `i` is
`(coords[i, j, 1], coords[i, j, 2])`, for `j in 1:n`, where point `n` repeats point 1. This is
the layout of `FixedWidthFloes.poly` (see `_floe_ring`), so any buffer with that layout, e.g.
one holding the overlap regions of two floes, can be used the same way. It's isbits on the
device, so kernels can create and pass it.
=#
struct _Ring{A<:AbstractArray{<:Any, 3}}
    coords::A
    i::Int
    n::Int
end

_npoint(r::_Ring) = r.n
@inline _getpoint(r::_Ring, j) = (r.coords[r.i, j, 1], r.coords[r.i, j, 2])

# Results of the point-in-polygon test, as GeometryOps' `point_in`, `point_on`, `point_out`
const _POINT_IN = Int8(1)
const _POINT_ON = Int8(0)
const _POINT_OUT = Int8(-1)

# Kinds of intersection between two edges, with the values of GeometryOps' `line_cross`,
# `line_hinge`, `line_over` and `line_out`
const _LINE_CROSS = Int8(1)
const _LINE_HINGE = Int8(2)
const _LINE_OVER = Int8(3)
const _LINE_OUT = Int8(4)

# Bounding box of the ring as (xmin, xmax, ymin, ymax), like `GI.extent`
function _ring_extent(r::_Ring)
    xmin, ymin = _getpoint(r, 1)
    xmax, ymax = xmin, ymin
    for j in 2:_npoint(r)
        x, y = _getpoint(r, j)
        xmin, xmax = min(xmin, x), max(xmax, x)
        ymin, ymax = min(ymin, y), max(ymax, y)
    end
    return xmin, xmax, ymin, ymax
end

# Whether the closed boxes (xmin, xmax, ymin, ymax) overlap or touch, like `Extents.intersects`
_extents_intersect(a, b) = a[1] <= b[2] && b[1] <= a[2] && a[3] <= b[4] && b[3] <= a[4]

#= Signed area of the ring, as `GO._signed_area(T, ring)`: the shoelace formula, closing the
ring from the last point to the first (which adds zero for a closed ring). =#
function _ring_signed_area(::Type{T}, r::_Ring) where T
    area = zero(T)
    n = _npoint(r)
    n == 0 && return area
    xfirst, yfirst = _getpoint(r, 1)
    x1, y1 = xfirst, yfirst
    for j in 2:n
        x2, y2 = _getpoint(r, j)
        area += x1 * y2 - y1 * x2
        x1, y1 = x2, y2
    end
    area += x1 * yfirst - y1 * xfirst
    return T(area / 2)
end

# Area of a polygon without holes, as `GO.area(poly, T)` (`area_poly`)
_ring_area(::Type{T}, r::_Ring) where T = abs(_ring_signed_area(T, r))

#= Centroid of a polygon without holes, as `GO.centroid(poly, T)` (`centroid_poly`). Like
GeometryOps, it computes the ring's centroid and then weights it by the area and divides by
the area again, as for a polygon with holes, which can change the last bit. =#
function _ring_centroid(::Type{T}, r::_Ring) where T
    xcentroid = zero(T)
    ycentroid = zero(T)
    area = zero(T)
    x1, y1 = _getpoint(r, 1)
    for j in 1:_npoint(r)
        x2, y2 = _getpoint(r, j)
        area_component = x1 * y2 - x2 * y1
        area += area_component
        xcentroid += (x1 + x2) * area_component
        ycentroid += (y1 + y2) * area_component
        x1, y1 = x2, y2
    end
    area /= 2
    xcentroid /= 6area
    ycentroid /= 6area
    area = abs(area)
    xcentroid *= area
    ycentroid *= area
    xcentroid /= area
    ycentroid /= area
    return xcentroid, ycentroid
end

#= Squared distance from (x0, y0) to the segment from (x1, y1) to (x2, y2), as
`GO._squared_distance_line`. The point and the segment can be of different types; the
result is converted to T at the end. =#
function _squared_distance_line(::Type{T}, x0, y0, x1, y1, x2, y2) where T
    xfirst, yfirst, xlast, ylast = x1 < x2 ? (x1, y1, x2, y2) : (x2, y2, x1, y1)
    v1, v2 = xlast - xfirst, ylast - yfirst
    w1, w2 = x0 - xfirst, y0 - yfirst
    c1 = w1 * v1 + w2 * v2
    c1 <= 0 && return T((xfirst - x0)^2 + (yfirst - y0)^2)
    c2 = v1 * v1 + v2 * v2
    c2 <= c1 && return T((xlast - x0)^2 + (ylast - y0)^2)
    b2 = c1 / c2
    return T((xfirst + (b2 * v1) - x0)^2 + (yfirst + (b2 * v2) - y0)^2)
end

_distance_line(::Type{T}, x0, y0, x1, y1, x2, y2) where T =
    sqrt(_squared_distance_line(T, x0, y0, x1, y1, x2, y2))

#= Distance from (x, y) to the ring's edges, as `GO._distance_curve(T, point, ring;
close_curve = true)`. =#
function _ring_distance(::Type{T}, x, y, r::_Ring) where T
    n = _npoint(r)
    first_last_equal = _getpoint(r, 1) == _getpoint(r, n)
    close_curve = first_last_equal
    n -= first_last_equal ? 1 : 0
    min_dist = typemax(T)
    x1, y1 = _getpoint(r, close_curve ? n : 1)
    for j in (close_curve ? 1 : 2):n
        x2, y2 = _getpoint(r, j)
        dist = _distance_line(T, x, y, x1, y1, x2, y2)
        min_dist = dist < min_dist ? dist : min_dist
        x1, y1 = x2, y2
    end
    return min_dist
end

#= Signed distance from (x, y) to a polygon without holes, as `GO.signed_distance(point,
poly, T)` (`dist_to_poly`): negative if the point is inside, i.e. within the polygon. =#
function _ring_signed_distance(::Type{T}, x, y, r::_Ring) where T
    dist = _ring_distance(T, x, y, r)
    return _point_within_ring(x, y, r) ? -dist : dist
end

#= Orientation of c relative to the line from a to b, without exact arithmetic, as
`GO.Predicates.orient(a, b, c; exact = False())`: 1 if c is to the left, -1 if it's to the
right, 0 if the cross product is approximately zero. Like GeometryOps, it converts the
points to Float64. =#
function _orient_inexact(a, b, c)
    ax, ay = Float64(a[1]) - Float64(c[1]), Float64(a[2]) - Float64(c[2])
    bx, by = Float64(b[1]) - Float64(c[1]), Float64(b[2]) - Float64(c[2])
    c_t1 = ax * by
    c_t2 = ay * bx
    return isapprox(c_t1, c_t2) ? 0.0 : sign(c_t1 - c_t2)
end

#= Orientation of c relative to the line from a to b, with exact arithmetic, as
`GO.Predicates.orient(a, b, c; exact = True())`: AdaptivePredicates' `orient2p` in Float64.
=#
_orient_exact(a, b, c) = AdaptivePredicates.orient2p(
    (Float64(a[1]), Float64(a[2])),
    (Float64(b[1]), Float64(b[2])),
    (Float64(c[1]), Float64(c[2])),
)

#= Whether (x, y) is inside, on or outside the ring, as
`GO._point_filled_curve_orientation(point, ring; exact = False())`: a ray-crossing test
that also detects points on an edge. =#
function _point_ring_orientation(x, y, r::_Ring)
    n = _npoint(r)
    n -= _getpoint(r, 1) == _getpoint(r, n) ? 1 : 0
    k = 0  # number of ray crossings
    xs, ys = _getpoint(r, n)
    for j in 1:n
        xe, ye = _getpoint(r, j)
        v1 = ys - y
        v2 = ye - y
        if !((v1 < 0 && v2 < 0) || (v1 > 0 && v2 > 0))
            u1, u2 = xs - x, xe - x
            f = _orient_inexact((xs, ys), (xe, ye), (x, y))
            if v2 > 0 && v1 ≤ 0
                f == 0 && return _POINT_ON
                f > 0 && (k += 1)
            elseif v1 > 0 && v2 ≤ 0
                f == 0 && return _POINT_ON
                f < 0 && (k += 1)
            elseif v2 == 0 && v1 < 0
                f == 0 && return _POINT_ON
            elseif v1 == 0 && v2 < 0
                f == 0 && return _POINT_ON
            elseif v1 == 0 && v2 == 0
                u2 ≤ 0 && u1 ≥ 0 && return _POINT_ON
                u1 ≤ 0 && u2 ≥ 0 && return _POINT_ON
            end
        end
        xs, ys = xe, ye
    end
    return iseven(k) ? _POINT_OUT : _POINT_IN
end

#= Orientation of (x, y) relative to a polygon without holes, as GeometryOps'
`_point_polygon_process`: a point outside the polygon's bounding box is outside, without
the ray-crossing test. =#
function _point_polygon_orientation(x, y, r::_Ring)
    xmin, xmax, ymin, ymax = _ring_extent(r)
    _extents_intersect((x, x, y, y), (xmin, xmax, ymin, ymax)) || return _POINT_OUT
    return _point_ring_orientation(x, y, r)
end

# Whether (x, y) is inside or on a polygon without holes, as `GO.coveredby` (`coveredby_poly`)
_point_coveredby_ring(x, y, r::_Ring) = _point_polygon_orientation(x, y, r) != _POINT_OUT

# Whether (x, y) is strictly inside a polygon without holes, as `GO.within`
_point_within_ring(x, y, r::_Ring) = _point_polygon_orientation(x, y, r) == _POINT_IN

# x / y clamped between ϵ and 1 - ϵ, as `GO._clamped_frac`
_clamped_frac(x::T, y::T, ϵ = zero(T)) where T = clamp(x / y, ϵ, one(T) - ϵ)

# Distance between two points, as `GO.distance(p1, p2, T)`
_point_distance(::Type{T}, p1, p2) where T = sqrt(T((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2))

# Whether p is in the closed box (xmin, xmax, ymin, ymax), as `GO._point_in_extent`
_point_in_extent(p, ext) = ext[1] ≤ p[1] ≤ ext[2] && ext[3] ≤ p[2] ≤ ext[4]

#=
Intersection of the edges (a1, a2) and (b1, b2), whose points are tuples of type T, as
`GO._intersection_point(T, (a1, a2), (b1, b2); exact = True())`. Returns the kind of
intersection (`_LINE_OUT`, `_LINE_HINGE`, `_LINE_CROSS` or `_LINE_OVER`) and two
intersection points with their fractions along each edge, `(point, (α, β))`. The second is
only used for `_LINE_OVER` (collinear edges that overlap).
=#
function _intersection_point(::Type{T}, a1, a2, b1, b2) where T
    line_orient = _LINE_OUT
    intr1 = ((zero(T), zero(T)), (zero(T), zero(T)))
    intr2 = intr1
    no_intr_result = (line_orient, intr1, intr2)
    a_ext = (minmax(a1[1], a2[1])..., minmax(a1[2], a2[2])...)
    b_ext = (minmax(b1[1], b2[1])..., minmax(b1[2], b2[2])...)
    _extents_intersect(a_ext, b_ext) || return no_intr_result
    a1_orient = _orient_exact(b1, b2, a1)
    a2_orient = _orient_exact(b1, b2, a2)
    a1_orient != 0 && a1_orient == a2_orient && return no_intr_result
    b1_orient = _orient_exact(a1, a2, b1)
    b2_orient = _orient_exact(a1, a2, b2)
    b1_orient != 0 && b1_orient == b2_orient && return no_intr_result
    if a1_orient == a2_orient == b1_orient == b2_orient == 0
        line_orient, intr1, intr2 =
            _find_collinear_intersection(T, a1, a2, b1, b2, a_ext, b_ext, no_intr_result)
    elseif a1_orient == 0 || a2_orient == 0 || b1_orient == 0 || b2_orient == 0
        line_orient = _LINE_HINGE
        intr1 = _find_hinge_intersection(T, a1, a2, b1, b2, a1_orient, a2_orient, b1_orient)
    else
        line_orient = _LINE_CROSS
        intr1 = _find_cross_intersection(T, a1, a2, b1, b2, a_ext, b_ext)
    end
    return line_orient, intr1, intr2
end

# Overlap of two collinear edges, as `GO._find_collinear_intersection`
function _find_collinear_intersection(::Type{T}, a1, a2, b1, b2, a_ext, b_ext, no_intr_result) where T
    line_orient, intr1, intr2 = no_intr_result
    a1_in_b = _point_in_extent(a1, b_ext)
    a2_in_b = _point_in_extent(a2, b_ext)
    b1_in_a = _point_in_extent(b1, a_ext)
    b2_in_a = _point_in_extent(b2, a_ext)
    a_dist, b_dist = _point_distance(T, a1, a2), _point_distance(T, b1, b2)
    if a1_in_b && a2_in_b
        line_orient = _LINE_OVER
        β1 = _clamped_frac(_point_distance(T, a1, b1), b_dist)
        β2 = _clamped_frac(_point_distance(T, a2, b1), b_dist)
        intr1 = (a1, (zero(T), β1))
        intr2 = (a2, (one(T), β2))
    elseif b1_in_a && b2_in_a
        line_orient = _LINE_OVER
        α1 = _clamped_frac(_point_distance(T, b1, a1), a_dist)
        α2 = _clamped_frac(_point_distance(T, b2, a1), a_dist)
        intr1 = (b1, (α1, zero(T)))
        intr2 = (b2, (α2, one(T)))
    elseif a1_in_b && b1_in_a
        if a1 == b1
            line_orient = _LINE_HINGE
            intr1 = (a1, (zero(T), zero(T)))
        else
            line_orient = _LINE_OVER
            intr1, intr2 = _set_ab_collinear_intrs(T, a1, b1, zero(T), zero(T), a1, b1, a_dist, b_dist)
        end
    elseif a1_in_b && b2_in_a
        if a1 == b2
            line_orient = _LINE_HINGE
            intr1 = (a1, (zero(T), one(T)))
        else
            line_orient = _LINE_OVER
            intr1, intr2 = _set_ab_collinear_intrs(T, a1, b2, zero(T), one(T), a1, b1, a_dist, b_dist)
        end
    elseif a2_in_b && b1_in_a
        if a2 == b1
            line_orient = _LINE_HINGE
            intr1 = (a2, (one(T), zero(T)))
        else
            line_orient = _LINE_OVER
            intr1, intr2 = _set_ab_collinear_intrs(T, a2, b1, one(T), zero(T), a1, b1, a_dist, b_dist)
        end
    elseif a2_in_b && b2_in_a
        if a2 == b2
            line_orient = _LINE_HINGE
            intr1 = (a2, (one(T), one(T)))
        else
            line_orient = _LINE_OVER
            intr1, intr2 = _set_ab_collinear_intrs(T, a2, b2, one(T), one(T), a1, b1, a_dist, b_dist)
        end
    end
    return line_orient, intr1, intr2
end

# As `GO._set_ab_collinear_intrs`
_set_ab_collinear_intrs(::Type{T}, a_pt, b_pt, a_pt_α, b_pt_β, a1, b1, a_dist, b_dist) where T =
    (
        (a_pt, (a_pt_α, _clamped_frac(_point_distance(T, a_pt, b1), b_dist))),
        (b_pt, (_clamped_frac(_point_distance(T, b_pt, a1), a_dist), b_pt_β)),
    )

# Edges that touch at one point but aren't collinear, as `GO._find_hinge_intersection`
function _find_hinge_intersection(::Type{T}, a1, a2, b1, b2, a1_orient, a2_orient, b1_orient) where T
    pt, α, β = if a1 == b1
        a1, zero(T), zero(T)
    elseif a1 == b2
        a1, zero(T), one(T)
    elseif a2 == b1
        a2, one(T), zero(T)
    elseif a2 == b2
        a2, one(T), one(T)
    elseif a1_orient == 0
        a1, zero(T), _clamped_frac(_point_distance(T, b1, a1), _point_distance(T, b1, b2), eps(T))
    elseif a2_orient == 0
        a2, one(T), _clamped_frac(_point_distance(T, b1, a2), _point_distance(T, b1, b2), eps(T))
    elseif b1_orient == 0
        b1, _clamped_frac(_point_distance(T, a1, b1), _point_distance(T, a1, a2), eps(T)), zero(T)
    else  # b2_orient == 0
        b2, _clamped_frac(_point_distance(T, a1, b2), _point_distance(T, a1, a2), eps(T)), one(T)
    end
    return pt, (α, β)
end

# Edges that cross at a point that isn't an end of either, as `GO._find_cross_intersection`
function _find_cross_intersection(::Type{T}, a1, a2, b1, b2, a_ext, b_ext) where T
    (a1x, a1y), (a2x, a2y) = a1, a2
    Δax, Δay = a2x - a1x, a2y - a1y
    (b1x, b1y), (b2x, b2y) = b1, b2
    Δbx, Δby = b2x - b1x, b2y - b1y
    Δbax = b1x - a1x
    Δbay = b1y - a1y
    a_cross_b = Δax * Δby - Δay * Δbx
    α = _clamped_frac(Δbax * Δby - Δbay * Δbx, a_cross_b, eps(T))
    β = _clamped_frac(Δbax * Δay - Δbay * Δax, a_cross_b, eps(T))
    x = if Δax == 0
        a1x
    elseif Δbx == 0
        b1x
    else
        (a1x + α * Δax + b1x + β * Δbx) / 2
    end
    y = if Δay == 0
        a1y
    elseif Δby == 0
        b1y
    else
        (a1y + α * Δay + b1y + β * Δby) / 2
    end
    pt = (x, y)
    if !_point_in_extent(pt, a_ext) || !_point_in_extent(pt, b_ext)
        pt, α, β = _nearest_endpoint(T, a1, a2, b1, b2)
    end
    return pt, (α, β)
end

#= End of either edge that is closest to the other edge, as `GO._nearest_endpoint`. Like
GeometryOps, the fraction for b2 is measured from a2, not a1. =#
function _nearest_endpoint(::Type{T}, a1, a2, b1, b2) where T
    a_dist = _point_distance(T, a1, a2)
    b_dist = _point_distance(T, b1, b2)
    min_pt, min_dist = a1, _distance_line(T, a1..., b1..., b2...)
    α, β = eps(T), _clamped_frac(_point_distance(T, min_pt, b1), b_dist, eps(T))
    dist = _distance_line(T, a2..., b1..., b2...)
    if dist < min_dist
        min_pt, min_dist = a2, dist
        α, β = one(T) - eps(T), _clamped_frac(_point_distance(T, min_pt, b1), b_dist, eps(T))
    end
    dist = _distance_line(T, b1..., a1..., a2...)
    if dist < min_dist
        min_pt, min_dist = b1, dist
        α, β = _clamped_frac(_point_distance(T, min_pt, a1), a_dist, eps(T)), eps(T)
    end
    dist = _distance_line(T, b2..., a1..., a2...)
    if dist < min_dist
        min_pt, min_dist = b2, dist
        α, β = _clamped_frac(_point_distance(T, min_pt, a2), a_dist, eps(T)), one(T) - eps(T)
    end
    return min_pt, α, β
end

#=
Insert point p into the sorted list of distinct points in `out[k, 1:m, :]`, keeping it sorted
and distinct, as `unique!(sort!(points))` would. Returns the new number of points, or -1 if
there is no room.
=#
function _insert_sorted_unique!(out, k, m, p)
    pos = m + 1
    for l in 1:m
        q = (out[k, l, 1], out[k, l, 2])
        isequal(q, p) && return m
        if isless(p, q)
            pos = l
            break
        end
    end
    m == size(out, 2) && return -1
    for l in m:-1:pos
        out[k, l + 1, 1] = out[k, l, 1]
        out[k, l + 1, 2] = out[k, l, 2]
    end
    out[k, pos, 1] = p[1]
    out[k, pos, 2] = p[2]
    return m + 1
end

#=
Points where the edges of two rings intersect, as `GO.intersection_points(p1, p2, T)`
(`get_intersection_points`, with T = Float64): sorted and distinct. They're written to
`out[k, 1:m, :]`, and the number of points m is returned, or -1 if `out` has too few columns
for them.
=#
function _ring_intersection_points!(out, k, ::Type{T}, r1::_Ring, r2::_Ring) where T
    m = 0
    _extents_intersect(_ring_extent(r1), _ring_extent(r2)) || return m
    for i in 1:(_npoint(r1) - 1)
        p1, p2 = _getpoint(r1, i), _getpoint(r1, i + 1)
        a1, a2 = (T(p1[1]), T(p1[2])), (T(p2[1]), T(p2[2]))
        for j in 1:(_npoint(r2) - 1)
            q1, q2 = _getpoint(r2, j), _getpoint(r2, j + 1)
            b1, b2 = (T(q1[1]), T(q1[2])), (T(q2[1]), T(q2[2]))
            line_orient, (pt1, _), (pt2, _) = _intersection_point(T, a1, a2, b1, b2)
            line_orient == _LINE_OUT && continue
            m = _insert_sorted_unique!(out, k, m, pt1)
            m < 0 && return m
            if line_orient == _LINE_OVER
                m = _insert_sorted_unique!(out, k, m, pt2)
                m < 0 && return m
            end
        end
    end
    return m
end

#=
Indices of the ring's points that are within `atol` of each of the points
`points[k, 1:m, :]`, in increasing order, as `which_vertices_match_points(points, region,
atol)`. For each point, it takes the closest of the ring's points, the first one if there is
a tie; an index can appear more than once. Like that function:
- the distances are square roots of GeometryOps distances, so a point matches if it's
  closer than `atol^2`;
- if the first and last points are equal, the last isn't used, so a single point matches
  nothing.
The indices are written to `out[k, 1:c]`, and their number c is returned, or -1 if `out` has
too few columns for them.
=#
function _which_vertices_match_points!(out, k, points, m, r::_Ring, atol = 1)
    FT = eltype(r.coords)
    c = 0
    m == 0 && return c  # which_vertices_match_points throws a BoundsError
    npoints = m
    if (points[k, 1, 1], points[k, 1, 2]) == (points[k, m, 1], points[k, m, 2])
        npoints -= 1
    end
    for l in 1:npoints
        p = (points[k, l, 1], points[k, l, 2])
        min_dist = FT(Inf)
        min_vert = 1
        for j in 1:_npoint(r)
            dist = sqrt(_point_distance(FT, _getpoint(r, j), p))
            if dist < min_dist
                min_dist = dist
                min_vert = j
            end
        end
        if min_dist < atol
            c == size(out, 2) && return -1
            # insertion sort
            pos = c + 1
            while pos > 1 && out[k, pos - 1] > min_vert
                out[k, pos] = out[k, pos - 1]
                pos -= 1
            end
            out[k, pos] = min_vert
            c += 1
        end
    end
    return c
end
