"""
    centroid(poly::ArrayPolygon, ::Type{T})

Centroid of the polygon as a tuple `(x, y)` of type `T`, as `GO.centroid(poly, T)`.
"""
function centroid(p::ArrayPolygon, ::Type{T}) where T
    xcentroid = zero(T)
    ycentroid = zero(T)
    area = zero(T)
    x1, y1 = getpoint(p, 1)
    for j in 1:npoint(p)
        x2, y2 = getpoint(p, j)
        area_component = x1 * y2 - x2 * y1
        area += area_component
        xcentroid += (x1 + x2) * area_component
        ycentroid += (y1 + y2) * area_component
        x1, y1 = x2, y2
    end
    area /= 2
    xcentroid /= 6area
    ycentroid /= 6area
    #= Like GeometryOps, which computes the ring's centroid and then weights it by the area
    and divides by the area again, as for a polygon with holes. This can change the last
    bit. =#
    area = abs(area)
    xcentroid *= area
    ycentroid *= area
    xcentroid /= area
    ycentroid /= area
    return xcentroid, ycentroid
end
