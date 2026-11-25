"""
    ArrayPolygon(coords, i, n)

A polygon without holes, stored in row `i` of a 3D array: point `j` is
`(coords[i, j, 1], coords[i, j, 2])`, for `j in 1:n`, and point `n` repeats point 1, so the
ring is closed. Rows can hold many polygons, e.g. one per floe, padded to the same width.

It doesn't copy `coords`, so functions that change the polygon (e.g.
[`transform!`](@ref GeometryOpsGPU.transform!)) change the array. It's isbits on a device,
so kernels can create one and pass it to other functions.
"""
struct ArrayPolygon{A<:AbstractArray{<:Any, 3}}
    coords::A
    i::Int
    n::Int
end

"""
    npoint(poly::ArrayPolygon)

Number of points of the polygon's ring, including the last one that repeats the first, as
`GI.npoint`.
"""
npoint(p::ArrayPolygon) = p.n

"""
    getpoint(poly::ArrayPolygon, j)

Point `j` of the polygon's ring as a tuple `(x, y)`, as `GI.getpoint(poly, j)`.
"""
@inline getpoint(p::ArrayPolygon, j) = (p.coords[p.i, j, 1], p.coords[p.i, j, 2])

# Set point j of the polygon's ring to (x, y)
@inline function _setpoint!(p::ArrayPolygon, j, (x, y))
    p.coords[p.i, j, 1] = x
    p.coords[p.i, j, 2] = y
    return
end
