"""
    area(poly::ArrayPolygon, ::Type{T})

Area of the polygon, in type `T`, as `GO.area(poly, T)`.
"""
area(p::ArrayPolygon, ::Type{T}) where T = abs(_signed_area(T, p))

#= Signed area of the polygon's ring, as `GO._signed_area(T, ring)`: the shoelace formula,
closing the ring from the last point to the first (which adds zero for a closed ring). =#
function _signed_area(::Type{T}, p::ArrayPolygon) where T
    area = zero(T)
    n = npoint(p)
    n == 0 && return area
    xfirst, yfirst = getpoint(p, 1)
    x1, y1 = xfirst, yfirst
    for j in 2:n
        x2, y2 = getpoint(p, j)
        area += x1 * y2 - y1 * x2
        x1, y1 = x2, y2
    end
    area += x1 * yfirst - y1 * xfirst
    return T(area / 2)
end
