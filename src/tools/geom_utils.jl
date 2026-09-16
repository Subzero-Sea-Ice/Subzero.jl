export make_polygon

#=
Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates
 are the exterior border of the floe and the wz coordinates, or any other
 following sets of coordinates, describe holes within the floe.
 This form is for easy conversion to polygons.
=#
const PolyVec{T} = Vector{Vector{Vector{T}}} where T<:Real

#=
Coordinates are vector of vector of points of the form:
[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]] where the xy coordinates form a
closed ring. PolyVec objects can be made out RingVec objects.
This form is for each conversion to LinearRings, which can also be made into Polygons.
=#
const RingVec{T} = R where {
    T<:Real,
    V<:AbstractArray{T},
    R <: AbstractArray{V},
}

# Define very specific type that GeometryOps returns so that it can be used to dispatch within Subzerp
const Polys{T} = GI.Polygon{false, false, Vector{GI.LinearRing{false, false, Vector{Tuple{T, T}}, Nothing, Nothing}}, Nothing, Nothing} where T
const MultiPolys{T} = GI.MultiPolygon{false, false, Vector{Polys{T}}, Nothing, Nothing} where T
const StaticQuadrilateral{FT} =  GI.Polygon{false,false, SA.SVector{1, GI.LinearRing{false, false, SA.SVector{5, Tuple{FT, FT}}, Nothing, Nothing}},Nothing,Nothing} where FT

# Convert polygons with points of type Float32/Float64 to type Float64/Float32
Base.convert(::Type{Polys{Float32}}, p::Polys{<:Real}) = GO.tuples(p, Float32)
Base.convert(::Type{Polys{Float64}}, p::Polys{<:Real}) = GO.tuples(p, Float64)

#=
Takes a RingVec object and make sure that the last element has the same first
element as last element and that other than these two elements there are no
duplicate, adjacent vertices. Also asserts that the ring as at least three
elements or else it cannot be made into a valid ring as it is a line segment. 
=#
function valid_ringvec!(ring)
    deleteat!(ring, findall(i->ring[i]==ring[i+1], collect(1:length(ring)-1)))
    if ring[1] != ring[end]
        push!(ring, deepcopy(ring[1]))
    end
    @assert length(ring) > 3 "Polgon needs at least 3 distinct points."
    return ring
end

#=
Takes a PolyVec object and make sure that the last element of each "ring"
(vector of vector of floats) has the same first element as last element and has
not duplicate adjacent elements. Also asserts that each "ring" as at least three
distinct elements or else it is not a valid ring, but rather a line segment. 
=#
function valid_polyvec!(coords)
    for ring in coords
        valid_ringvec!(ring)
    end
    return coords
end

# Wrappers for calling GeometryOps (GO) functions within the code!!! If you wanted to use a
# different library or dispatch for a specific type of polygon (i.e. disks) then you would need
# to re-write most of these (unless you don't want the specific functionality they offer, like fracturing)

# find the intersection of two polygons and return as a list of polygons
intersect_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.intersection(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing)
# find the difference of two polygons and return as a list of polygons
diff_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.difference(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing) 
# find the union of two polygons and return as a list of polygons
union_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.union(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing)
# simplify an existing polygon to have less vertices
simplify_poly(p, tol) = GO.simplify(p; tol = tol)
# area of a polygon
area_poly(p, ::Type{FT}) where FT = GO.area(p, FT)
# centroid of a polygon
centroid_poly(p, ::Type{FT}) where FT = GO.centroid(p, FT)
# signed distance from a point to polygon
dist_to_poly(point, poly, ::Type{FT}) where FT = GO.signed_distance(point, poly, FT)
# boolean if point is covered by a polygon
coveredby_poly(point, poly) = GO.coveredby(point, poly)
# return point as a tuple with elements of type FT
get_tuple_point(point, ::Type{FT}) where FT = GO._tuple_point(point, FT)
# return polygon with points represented as tuples of type FT
get_tuple_poly(poly, ::Type{FT}) where FT = GO.tuples(poly, FT)
# get internal angles of polygon
angles_poly(poly, ::Type{FT}) where FT = GO.angles(poly, FT)
# check intersections between polys
check_intersects(poly1, poly2) = GO.intersects(poly1, poly2)
# get intersection points between polys
get_intersection_points(poly1, poly2) = GO.intersection_points(poly1, poly2)
# cut polygon by line through it and return polys of type FT
cut_poly_by_line(poly, line, ::Type{FT}) where FT = GO.cut(poly, line, FT)

# translate polygon coordinates by Δx, Δy
function _translate_poly(::Type{FT}, p, Δx, Δy) where FT
    t = CoordinateTransformations.Translation(Δx, Δy)
    return get_tuple_poly(GO.transform(t, p), FT)
end
# translate polygon coordinates by Δx, Δy and rotate polygon coordiantes by Δα
function _move_poly(poly::Subzero.Polys{FT}, Δx, Δy, Δα, cx = 0.0, cy = 0.0) where FT
    rot = CoordinateTransformations.LinearMap(Rotations.Angle2d(Δα))
    cent_rot = CoordinateTransformations.recenter(rot, (cx, cy))
    trans = CoordinateTransformations.Translation(Δx, Δy)
    # TODO: can remove the tuples call after GO SVPoint PR
    return get_tuple_poly(GO.transform(trans ∘ cent_rot, poly), FT)
end

# create polygon from a PolyVec, tuple coordiantes, or a linear ring
make_polygon(coords, ::Type{FT} = Float64) where FT = GI.Polygon(get_tuple_poly(coords, FT))
# make_polygon(tuple_coords, ::Type{FT}) = GI.Polygon(get_tuple_poly(coords, FT))
make_polygon(ring::GI.LinearRing, ::Type{FT} = Float64) where FT = GI.Polygon([get_tuple_poly(ring, FT)])
# create a multipolygon from a vector of PolyVecs, tuple coordiantes, a vector of polygons, or a
# vector of StaticQuadrilaterals (used for bounding boxes)
make_multipolygon(coords::Vector{<:PolyVec}, ::Type{FT} = Float64) where FT = GI.MultiPolygon(get_tuple_poly(coords, FT))
# make_multipolygon(tuple_coords, ::) = GI.MultiPolygon(tuple_coords)
make_multipolygon(polys::Vector{<:GI.Polygon}, ::Type{FT} = Float64) where FT = GI.MultiPolygon(get_tuple_poly.(polys, FT))
function make_multipolygon(polys::Vector{<:StaticQuadrilateral{FT}}) where FT
    new_polys = Vector{Polys{FT}}(undef, length(polys))
    for (i, poly) in enumerate(polys)
        new_polys[i] = make_polygon([[p for p in GI.getpoint(poly)]], FT)
    end
    return GI.MultiPolygon(new_polys)
end

# create a bounding box polygon (a rectangle!) that is used to create domain boundaries
function _make_bounding_box_polygon(::Type{FT}, xmin, xmax, ymin, ymax) where FT
    points = ((xmin, ymin),  (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin))
    ring = GI.LinearRing(SA.SVector{5, Tuple{FT, FT}}(points))
    return  GI.Polygon(SA.SVector(ring))
end

# Determine if polygon has one or more holes
function hashole(poly::Polys)
    return GI.nhole(poly) > 0
end 
# Remove any existing holes from polygon
function rmholes!(poly::Polys)
    deleteat!(poly.geom, 2:GI.nring(poly))
    return
end

# Find the length of the maximum radius of a given polygon
function _calc_max_radius(poly, cent, ::Type{T}) where T
    max_rad_sqrd = zero(T)
    Δx, Δy = get_tuple_point(cent, T)
    for pt in GI.getpoint(GI.getexterior(poly))
        x, y = get_tuple_point(pt, T)
        x, y = x - Δx, y - Δy
        rad_sqrd = x^2 + y^2
        if rad_sqrd > max_rad_sqrd
            max_rad_sqrd = rad_sqrd
        end
    end
    return sqrt(max_rad_sqrd)
end

#=
Find which vertices in in a polygon match the user-provided points.
Provided points can be a vector of tuples or a vector of vectors.
Provide sorted indices which tell which of the polyon region's vertices
are within atol of any of the provided points.

Note: If last coordinate is a repeat of first coordinate, last coordinate
    index is NOT recorded.
=#
function which_vertices_match_points(points, region::Polys{FT}, atol = 1) where FT
    idxs = Vector{Int}()
    npoints = length(points)
    if points[1] == points[end]
        npoints -= 1
    end
    for i in 1:npoints  # find which vertex matches point
        min_dist = FT(Inf)
        min_vert = 1
        for (j, pt) in enumerate(GI.getpoint(GI.getexterior(region)))
            dist = sqrt(GO.distance(pt, points[i], FT))
            if dist < min_dist
                min_dist = dist
                min_vert = j
            end
        end
        if min_dist < atol
            push!(idxs, min_vert)
        end
    end
    return sort!(idxs)
end

#=
Generate voronoi coords within a bounding box defined by its lower left corner
and its height and width. Attempt to generate `npieces` cells within the box.
Inputs:
    desired_points  <Int> desired number of voronoi cells
    Δx              <AbstractFloat> width of bounding box
    Δy              <AbstractFloat> height of bounding box
    xmin            <AbstractFloat> minimum x-value of bounding box
    ymin            <AbstractFloat> minimum y-value of bounding box 
    domain_coords   <Vector{PolyVec{AbstractFloat}}> multipolygon that will
                        eventually be filled with/intersected with the voronoi
                        cells - such as topography
    rng             <RNG> random number generator to generate voronoi cells
    min_to_warn     <Int> minimum number of points to warn if not generated to
                        seed voronoi
    max_tries       <Int> number of tires to generate desired number of points
                        within domain_coords to seed voronoi cell creation
Outputs:
    coords  <Vector{PolyVec{Float}}> vector of polygon coordinates generated by
        voronoi tesselation. These polygons all fall within the space defined by
        the domain_coords. If less polygons than min_to_warn are generated, the
        user will be warned. 
=#
function _generate_voronoi_coords(::Type{FT}, desired_points::Int, Δx, Δy, xmin, ymin,
    domain_poly, rng, min_to_warn::Int; max_tries::Int = 10,
) where {FT <: AbstractFloat}
    xpoints = Vector{Float64}()
    ypoints = Vector{Float64}()
    area_frac = area_poly(domain_poly, FT) / (Δx * Δy)
    # Increase the number of points based on availible percent of bounding box
    npoints = ceil(Int, desired_points / area_frac)
    current_points = 0
    tries = 0
    while current_points < desired_points && tries <= max_tries
        x = (rand(rng, npoints) .* Δx) .+ xmin
        y = (rand(rng, npoints) .* Δy) .+ ymin
        # Check which of the scaled and translated points are within the domain coords
        in_idx = [coveredby_poly((x[i], y[i]), domain_poly) for i in eachindex(x)]
        current_points += sum(in_idx)
        tries += 1
        append!(xpoints, x[in_idx])
        append!(ypoints, y[in_idx])
    end
    # If we generated too many cells, remove extra
    if current_points > desired_points
        xpoints = xpoints[1:desired_points]
        ypoints = ypoints[1:desired_points]
        current_points = desired_points
    end
    # Warn if we didn't generate enough cells
    if current_points < min_to_warn
        @warn "Only $current_points floes were able to be generated in \
            $max_tries tries during voronoi tesselation."
    end
    # Make voronoi cells into floes
    polys = if current_points > 2
        xmax, ymax = xmin + Δx, ymin + Δy
        clip_points = GO.tuples(((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)), FT)
        clip_vertices = (1, 2, 3, 4, 1)
        clip_polygon = (clip_points, clip_vertices)
        try
            GO.voronoi(tuple.(xpoints, ypoints), FT; clip_polygon)
        catch
            @info "Voronoi failed - floe to fracture likely really small." tstep = tstep
            Polys{FT}[]
        end
    else
        Polys{FT}[]
    end
    return polys
end