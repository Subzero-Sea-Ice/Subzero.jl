"""
    distance(p1, p2, ::Type{T})

Distance between the points `p1` and `p2`, tuples `(x, y)` of type `T`, as
`GO.distance(p1, p2, T)`.
"""
distance(p1::NTuple{2, T}, p2::NTuple{2, T}, ::Type{T}) where T = _euclid_distance(T, p1, p2)

"""
    signed_distance(point, poly::ArrayPolygon, ::Type{T})

Distance from the point, a tuple `(x, y)`, to the polygon's edges, as
`GO.signed_distance(point, poly, T)`: negative if the point is within the polygon. The point
and `T` have the polygon's float type.
"""
function signed_distance(point::NTuple{2, T}, p::ArrayPolygon{<:AbstractArray{T}},
        ::Type{T}) where T
    dist = _distance_curve(T, point, p)
    return within(point, p) ? -dist : dist
end

# Distance between two points, as `GO._euclid_distance(T, p1, p2)`
_euclid_distance(::Type{T}, p1::NTuple{2, T}, p2::NTuple{2, T}) where T =
    sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)

# Squared distance from point p0 to the segment from p1 to p2, as `GO._squared_distance_line`
function _squared_distance_line(::Type{T}, p0::NTuple{2, T}, p1::NTuple{2, T},
        p2::NTuple{2, T}) where T
    x0, y0 = p0
    x1, y1 = p1
    x2, y2 = p2
    xfirst, yfirst, xlast, ylast = x1 < x2 ? (x1, y1, x2, y2) : (x2, y2, x1, y1)
    v1, v2 = xlast - xfirst, ylast - yfirst
    w1, w2 = x0 - xfirst, y0 - yfirst
    c1 = w1 * v1 + w2 * v2
    c1 <= 0 && return (xfirst - x0)^2 + (yfirst - y0)^2
    c2 = v1 * v1 + v2 * v2
    c2 <= c1 && return (xlast - x0)^2 + (ylast - y0)^2
    b2 = c1 / c2
    return (xfirst + (b2 * v1) - x0)^2 + (yfirst + (b2 * v2) - y0)^2
end

# Distance from point p0 to the segment from p1 to p2, as `GO._distance_line`
_distance_line(::Type{T}, p0, p1, p2) where T = sqrt(_squared_distance_line(T, p0, p1, p2))

#= Distance from the point to the polygon's edges, as `GO._distance_curve(T, point, ring;
close_curve = true)`. =#
function _distance_curve(::Type{T}, point::NTuple{2, T},
        p::ArrayPolygon{<:AbstractArray{T}}) where T
    n = npoint(p)
    first_last_equal = getpoint(p, 1) == getpoint(p, n)
    close_curve = first_last_equal
    n -= first_last_equal ? 1 : 0
    min_dist = typemax(T)
    p1 = getpoint(p, close_curve ? n : 1)
    for j in (close_curve ? 1 : 2):n
        p2 = getpoint(p, j)
        dist = _distance_line(T, point, p1, p2)
        min_dist = dist < min_dist ? dist : min_dist
        p1 = p2
    end
    return min_dist
end
