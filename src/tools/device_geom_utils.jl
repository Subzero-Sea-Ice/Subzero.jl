#=
Methods of the geometry functions in `geom_utils.jl` for polygons stored in arrays
(`GOG.ArrayPolygon`s, e.g. `_floe_poly(floes, i)`), so that they can be called in GPU
kernels. They call `GeometryOpsGPU`, as the methods in `geom_utils.jl` call GeometryOps.
Functions that write into a buffer instead of returning a vector end in `!`.
=#

area_poly(p::GOG.ArrayPolygon, ::Type{T}) where T = GOG.area(p, T)

centroid_poly(p::GOG.ArrayPolygon, ::Type{T}) where T = GOG.centroid(p, T)

dist_to_poly(point, p::GOG.ArrayPolygon, ::Type{T}) where T = GOG.signed_distance(point, p, T)

coveredby_poly(point, p::GOG.ArrayPolygon) = GOG.coveredby(point, p)

#=
Points where the edges of two polygons intersect, as `get_intersection_points(poly1, poly2)`:
sorted and distinct. They're written to `out[k, 1:m, :]`, and the number of points m is
returned, or -1 if `out` has too few columns for them.
=#
get_intersection_points!(out, k, poly1::GOG.ArrayPolygon, poly2::GOG.ArrayPolygon,
    ::Type{T} = Float64) where T = GOG.intersection_points!(out, k, poly1, poly2, T)

#=
Indices of the region's points that are within `atol` of each of the points
`points[k, 1:m, :]`, in increasing order, as `which_vertices_match_points(points, region,
atol)`. For each point, it takes the closest of the region's points, the first one if there
is a tie; an index can appear more than once. Like that function:
- the distances are square roots of distances, so a point matches if it's closer than
  `atol^2`;
- if the first and last points are equal, the last isn't used, so a single point matches
  nothing.
The indices are written to `out[k, 1:c]`, and their number c is returned, or -1 if `out` has
too few columns for them.
=#
function which_vertices_match_points!(out, k, points, m, region::GOG.ArrayPolygon, atol = 1)
    FT = eltype(region.coords)
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
        for j in 1:GOG.npoint(region)
            dist = sqrt(GOG.distance(GOG.getpoint(region, j), p, FT))
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
