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
sorted and distinct, in the polygons' float type, which `out` and `T` must have. They're
written to `out[k, 1:m, :]`, and the number of points m is returned, or -1 if `out` has too
few columns for them.
=#
get_intersection_points!(out, k, poly1::GOG.ArrayPolygon, poly2::GOG.ArrayPolygon,
    ::Type{T} = eltype(poly1.coords)) where T = GOG.intersection_points!(out, k, poly1, poly2, T)

#=
The regions where two polygons overlap, as `intersect_polys(poly_a, poly_b, FT)`. They're
written to `out` (a `GOG.PolyRegions`), using the node lists `lists` (a `GOG.ClipLists`), and
their number is returned, or `GOG.CLIP_TOO_SMALL` or `GOG.CLIP_FAILED` (see
`GOG.intersection!`). It clips in the float type of `out` and `lists`: use the floes' type,
so that Float32 floes are clipped with Float32 maths. `intersect_polys` always clips in
Float64, also for Float32 floes, because `GO.intersection` doesn't pass the float type on to
the clipping, so with Float32 the regions differ from it by round-off.
=#
intersect_polys!(out::GOG.PolyRegions, lists::GOG.ClipLists, poly_a::GOG.ArrayPolygon,
    poly_b::GOG.ArrayPolygon) = GOG.intersection!(out, lists, poly_a, poly_b)

#=
Indices of the region's points that are within `atol` of each of the points
`points[k, 1:m, :]`, in increasing order, as `which_vertices_match_points(points, region,
atol)`. For each point, it takes the closest of the region's points, the first one if there
is a tie; an index can appear more than once. If the first and last points are equal, the
last isn't used, so a single point matches nothing.
The indices are written to `out[k, 1:c]`, and their number c is returned, or -1 if `out` has
too few columns for them. `points` must have the region's float type.
=#
function which_vertices_match_points!(out, k, points, m, region::GOG.ArrayPolygon, atol = 1)
    FT = eltype(region.coords)
    c = 0
    m == 0 && return c  # with no points, nothing matches
    npoints = m
    if (points[k, 1, 1], points[k, 1, 2]) == (points[k, m, 1], points[k, m, 2])
        npoints -= 1
    end
    for l in 1:npoints
        p = (points[k, l, 1], points[k, l, 2])
        min_dist = FT(Inf)
        min_vert = 1
        for j in 1:GOG.npoint(region)
            #= TODO: this is the square root of the distance, so a point matches if it's
            closer than `atol^2`. That's the same for the default `atol = 1`, but probably not
            what was meant. Kept, as in which_vertices_match_points. =#
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
