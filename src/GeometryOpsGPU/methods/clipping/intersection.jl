"""
    LineOrientation

How two edges intersect: `line_cross` (at a point that isn't an end of either),
`line_hinge` (at an end), `line_over` (collinear and overlapping), or `line_out` (not at
all), as GeometryOps' `LineOrientation`.
"""
@enum LineOrientation line_cross=1 line_hinge=2 line_over=3 line_out=4

"""
    intersection_points!(out, k, poly1::ArrayPolygon, poly2::ArrayPolygon, ::Type{T})

Points where the edges of the two polygons intersect, sorted and distinct, as
`GO.intersection_points(poly1, poly2, T)`. Unlike that, they are written into
`out[k, 1:m, :]` (point `l` is `(out[k, l, 1], out[k, l, 2])`), and the number of points
`m` is returned, or -1 if `out` has too few columns for them. `out` and `T` have the
polygons' float type.
"""
function intersection_points!(out::AbstractArray{T, 3}, k,
        poly1::ArrayPolygon{<:AbstractArray{T}}, poly2::ArrayPolygon{<:AbstractArray{T}},
        ::Type{T}) where T
    m = 0
    Extents.intersects(extent(poly1), extent(poly2)) || return m
    for i in 1:(npoint(poly1) - 1)
        a1, a2 = getpoint(poly1, i), getpoint(poly1, i + 1)
        for j in 1:(npoint(poly2) - 1)
            b1, b2 = getpoint(poly2, j), getpoint(poly2, j + 1)
            line_orient, (pt1, _), (pt2, _) =
                _intersection_point(T, (a1, a2), (b1, b2); exact = True())
            line_orient == line_out && continue
            m = _insert_sorted_unique!(out, k, m, pt1)
            m < 0 && return m
            if line_orient == line_over
                m = _insert_sorted_unique!(out, k, m, pt2)
                m < 0 && return m
            end
        end
    end
    return m
end

#=
Insert point p into the sorted list of distinct points in `out[k, 1:m, :]`, keeping it sorted
and distinct, so that the result is the same as GeometryOps' `unique!(sort!(points))`.
Returns the new number of points, or -1 if there is no room.
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
Intersection of the edges (a1, a2) and (b1, b2), whose points are tuples of type T, as
`GO._intersection_point(Planar(), T, (a1, a2), (b1, b2); exact)`. Returns the kind of
intersection (a `LineOrientation`) and two intersection points with their fractions along
each edge, `(point, (α, β))`. The second is only used for `line_over` (collinear edges that
overlap). `exact` (`True()` or `False()`) is passed on to the orientation tests: the
clipping and `intersection_points!` use `True()`, `disjoint` uses `False()`.
=#
function _intersection_point(::Type{T}, (a1, a2)::NTuple{2, NTuple{2, T}},
        (b1, b2)::NTuple{2, NTuple{2, T}}; exact) where T
    line_orient = line_out
    intr1 = ((zero(T), zero(T)), (zero(T), zero(T)))
    intr2 = intr1
    no_intr_result = (line_orient, intr1, intr2)
    a_ext = Extents.Extent(X = minmax(a1[1], a2[1]), Y = minmax(a1[2], a2[2]))
    b_ext = Extents.Extent(X = minmax(b1[1], b2[1]), Y = minmax(b1[2], b2[2]))
    Extents.intersects(a_ext, b_ext) || return no_intr_result
    a1_orient = Predicates.orient(b1, b2, a1; exact)
    a2_orient = Predicates.orient(b1, b2, a2; exact)
    a1_orient != 0 && a1_orient == a2_orient && return no_intr_result
    b1_orient = Predicates.orient(a1, a2, b1; exact)
    b2_orient = Predicates.orient(a1, a2, b2; exact)
    b1_orient != 0 && b1_orient == b2_orient && return no_intr_result
    if a1_orient == a2_orient == b1_orient == b2_orient == 0
        line_orient, intr1, intr2 =
            _find_collinear_intersection(T, a1, a2, b1, b2, a_ext, b_ext, no_intr_result)
    elseif a1_orient == 0 || a2_orient == 0 || b1_orient == 0 || b2_orient == 0
        line_orient = line_hinge
        intr1 = _find_hinge_intersection(T, a1, a2, b1, b2, a1_orient, a2_orient, b1_orient)
    else
        line_orient = line_cross
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
    a_dist, b_dist = _euclid_distance(T, a1, a2), _euclid_distance(T, b1, b2)
    if a1_in_b && a2_in_b
        line_orient = line_over
        β1 = _clamped_frac(_euclid_distance(T, a1, b1), b_dist)
        β2 = _clamped_frac(_euclid_distance(T, a2, b1), b_dist)
        intr1 = (a1, (zero(T), β1))
        intr2 = (a2, (one(T), β2))
    elseif b1_in_a && b2_in_a
        line_orient = line_over
        α1 = _clamped_frac(_euclid_distance(T, b1, a1), a_dist)
        α2 = _clamped_frac(_euclid_distance(T, b2, a1), a_dist)
        intr1 = (b1, (α1, zero(T)))
        intr2 = (b2, (α2, one(T)))
    elseif a1_in_b && b1_in_a
        if a1 == b1
            line_orient = line_hinge
            intr1 = (a1, (zero(T), zero(T)))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(T, a1, b1, zero(T), zero(T), a1, b1, a_dist, b_dist)
        end
    elseif a1_in_b && b2_in_a
        if a1 == b2
            line_orient = line_hinge
            intr1 = (a1, (zero(T), one(T)))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(T, a1, b2, zero(T), one(T), a1, b1, a_dist, b_dist)
        end
    elseif a2_in_b && b1_in_a
        if a2 == b1
            line_orient = line_hinge
            intr1 = (a2, (one(T), zero(T)))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(T, a2, b1, one(T), zero(T), a1, b1, a_dist, b_dist)
        end
    elseif a2_in_b && b2_in_a
        if a2 == b2
            line_orient = line_hinge
            intr1 = (a2, (one(T), one(T)))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(T, a2, b2, one(T), one(T), a1, b1, a_dist, b_dist)
        end
    end
    return line_orient, intr1, intr2
end

# As `GO._set_ab_collinear_intrs`
_set_ab_collinear_intrs(::Type{T}, a_pt, b_pt, a_pt_α, b_pt_β, a1, b1, a_dist, b_dist) where T =
    (
        (a_pt, (a_pt_α, _clamped_frac(_euclid_distance(T, a_pt, b1), b_dist))),
        (b_pt, (_clamped_frac(_euclid_distance(T, b_pt, a1), a_dist), b_pt_β)),
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
        a1, zero(T), _clamped_frac(_euclid_distance(T, b1, a1), _euclid_distance(T, b1, b2), eps(T))
    elseif a2_orient == 0
        a2, one(T), _clamped_frac(_euclid_distance(T, b1, a2), _euclid_distance(T, b1, b2), eps(T))
    elseif b1_orient == 0
        b1, _clamped_frac(_euclid_distance(T, a1, b1), _euclid_distance(T, a1, a2), eps(T)), zero(T)
    else  # b2_orient == 0
        b2, _clamped_frac(_euclid_distance(T, a1, b2), _euclid_distance(T, a1, a2), eps(T)), one(T)
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
    a_dist = _euclid_distance(T, a1, a2)
    b_dist = _euclid_distance(T, b1, b2)
    min_pt, min_dist = a1, _distance_line(T, a1, b1, b2)
    α, β = eps(T), _clamped_frac(_euclid_distance(T, min_pt, b1), b_dist, eps(T))
    dist = _distance_line(T, a2, b1, b2)
    if dist < min_dist
        min_pt, min_dist = a2, dist
        α, β = one(T) - eps(T), _clamped_frac(_euclid_distance(T, min_pt, b1), b_dist, eps(T))
    end
    dist = _distance_line(T, b1, a1, a2)
    if dist < min_dist
        min_pt, min_dist = b1, dist
        α, β = _clamped_frac(_euclid_distance(T, min_pt, a1), a_dist, eps(T)), eps(T)
    end
    dist = _distance_line(T, b2, a1, a2)
    if dist < min_dist
        min_pt, min_dist = b2, dist
        α, β = _clamped_frac(_euclid_distance(T, min_pt, a2), a_dist, eps(T)), one(T) - eps(T)
    end
    return min_pt, α, β
end

# x / y clamped between ϵ and 1 - ϵ, as `GO._clamped_frac`
_clamped_frac(x::T, y::T, ϵ = zero(T)) where T = clamp(x / y, ϵ, one(T) - ϵ)

"""
    intersection!(out::PolyRegions, lists::ClipLists, poly_a::ArrayPolygon, poly_b::ArrayPolygon)

The regions where the polygons overlap, as `GO.intersection(poly_a, poly_b, T; target =
GI.PolygonTrait(), fix_multipoly = nothing)` with Foster-Hormann clipping, where `T` is the
polygons' float type, which `out` and `lists` must have too. Unlike that, the regions are written
into `out` (see [`PolyRegions`](@ref GeometryOpsGPU.PolyRegions)), using `lists` as scratch
memory, and their number is returned: 0 if the polygons don't overlap or only touch.
Returns `CLIP_TOO_SMALL` if `lists` or `out` is too small, and `CLIP_FAILED` where
GeometryOps throws an error.
"""
function intersection!(out::PolyRegions{<:AbstractArray{T, 3}},
        lists::ClipLists{<:_PolyNodes{<:AbstractArray{T, 3}}},
        poly_a::ArrayPolygon{<:AbstractArray{T}}, poly_b::ArrayPolygon{<:AbstractArray{T}}) where T
    a, b, a_idx = lists.a, lists.b, lists.a_idx
    na, n_idx, _ = _build_a_list!(a, a_idx, T, poly_a, poly_b)
    na < 0 && return CLIP_TOO_SMALL
    na < 2 && return CLIP_FAILED  # GeometryOps can't index a_list[end - 1]
    nb = _build_b_list!(b, a, a_idx, n_idx, T, poly_b)
    nb < 0 && return CLIP_TOO_SMALL
    _classify_crossing!(a, na, b, nb) || return CLIP_FAILED
    _flag_ent_exit!(a, na, poly_b) || return CLIP_FAILED
    _flag_ent_exit!(b, nb, poly_a) || return CLIP_FAILED
    n_cross = _index_crossing_intrs!(a, b, a_idx, n_idx)
    n_polys = _trace_polynodes!(out, a, na, b, nb, a_idx, n_cross)
    n_polys < 0 && return n_polys
    if n_polys == 0  # the boundaries don't cross: one polygon may be inside the other
        a_in_b, b_in_a = _find_non_cross_orientation(a, na, b, nb, poly_a, poly_b)
        if a_in_b || b_in_a
            out.max_polys < 1 && return CLIP_TOO_SMALL
            _copy_poly_to_region!(out, 1, a_in_b ? poly_a : poly_b, T) ||
                return CLIP_TOO_SMALL
            n_polys = 1
        end
    end
    return _remove_collinear_points!(out, n_polys)
end
