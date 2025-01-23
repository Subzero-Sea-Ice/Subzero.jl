
"""
    valid_ringvec(coords::RingVec{FT})

Takes a RingVec object and make sure that the last element has the same first
element as last element and that other than these two elements there are no
duplicate, adjacent vertices. Also asserts that the ring as at least three
elements or else it cannot be made into a valid ring as it is a line segment. 
"""
function valid_ringvec!(ring)
    deleteat!(ring, findall(i->ring[i]==ring[i+1], collect(1:length(ring)-1)))
    if ring[1] != ring[end]
        push!(ring, deepcopy(ring[1]))
    end
    @assert length(ring) > 3 "Polgon needs at least 3 distinct points."
    return ring
end

"""
    valid_polyvec(coords::PolyVec{FT})

Takes a PolyVec object and make sure that the last element of each "ring"
(vector of vector of floats) has the same first element as last element and has
not duplicate adjacent elements. Also asserts that each "ring" as at least three
distinct elements or else it is not a valid ring, but rather a line segment. 
"""
function valid_polyvec!(coords)
    for ring in coords
        valid_ringvec!(ring)
    end
    return coords
end

"""
    find_poly_coords(poly)

Syntactic sugar for to find a polygon's coordinates
Input:
    poly    <Polygon>
Output:
    <PolyVec> representing the floe's coordinates xy plane
"""
find_poly_coords(poly) = GI.coordinates(poly)

"""
    intersect_polys(p1, p2)

Intersect two geometries and return a list of polygons resulting.
Inputs:
    p1  <AbstractGeometry>
    p2  <AbstractGeometry>
Output:
    Vector of Polygons
"""
intersect_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.intersection(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing)
diff_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.difference(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing) 
union_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.union(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing)
simplify_poly(p, tol) = GO.simplify(p; tol = tol)

function _translate_poly(::Type{FT}, p, Δx, Δy) where FT
    t = CoordinateTransformations.Translation(Δx, Δy)
    # TODO: can remove the tuples call after GO SVPoint PR
    return GO.tuples(GO.transform(t, p), FT)
end

function _translate_floe!(::Type{FT}, floe, Δx, Δy) where FT
    translate!(floe.coords, Δx, Δy)
    floe.centroid[1] += Δx
    floe.centroid[2] += Δy
    floe.poly = _translate_poly(FT, floe.poly, Δx, Δy)
    return
end

function _move_poly(::Type{FT}, poly, Δx, Δy, Δα, cx = zero(FT), cy = zero(FT)) where FT
    rot = CoordinateTransformations.LinearMap(Rotations.Angle2d(Δα))
    cent_rot = CoordinateTransformations.recenter(rot, (cx, cy))
    trans = CoordinateTransformations.Translation(Δx, Δy)
    # TODO: can remove the tuples call after GO SVPoint PR
    return GO.tuples(GO.transform(trans ∘ cent_rot, poly), FT)::Polys{FT}
end

function _move_floe!(::Type{FT}, floe, Δx, Δy, Δα) where FT
    cx, cy = floe.centroid
    # Move coordinates and centroid
    translate!(floe.coords, -cx, -cy)
    rotate_radians!(floe.coords, Δα)
    floe.centroid[1] += Δx
    floe.centroid[2] += Δy
    translate!(floe.coords, cx + Δx, cy + Δy)
    # Move Polygon
    floe.poly = _move_poly(FT, floe.poly, Δx, Δy, Δα, cx, cy)::Polys{FT}
    return 
end

make_polygon(coords::PolyVec) = GI.Polygon(GO.tuples(coords))
make_polygon(tuple_coords) = GI.Polygon(tuple_coords)
make_polygon(ring::GI.LinearRing) = GI.Polygon([ring])
make_multipolygon(coords::Vector{<:PolyVec}) = GI.MultiPolygon(GO.tuples(coords))
make_multipolygon(tuple_coords) = GI.MultiPolygon(tuple_coords)
make_multipolygon(polys::Vector{<:GI.Polygon}) = GI.MultiPolygon(polys)

get_floe(floes::StructArray, i::Int) = LazyRow(floes, i)

function _make_bounding_box_polygon(::Type{FT}, xmin, xmax, ymin, ymax) where FT
    points = ((xmin, ymin),  (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin))
    ring = GI.LinearRing(SA.SVector{5, Tuple{FT, FT}}(points))
    return  GI.Polygon(SA.SVector(ring))
end

"""
    deepcopy_floe(floe::LazyRow{Floe{FT}})

Deepcopy of a floe by creating a new floe and copying all fields.
Inputs:
    floe    <Floe>
Outputs:
    New floe with floes that are equal in value. Any vector fields are copies so
    they share values, but not referance.
"""
function deepcopy_floe(floe::LazyRow{Floe{FT}}) where {FT}
    poly = GO.tuples(floe.poly, FT)
    f = Floe{FT}(
        poly = poly,
        centroid = copy(floe.centroid),
        coords = find_poly_coords(poly),
        height = floe.height,
        area = floe.area,
        mass = floe.mass,
        rmax = floe.rmax,
        moment = floe.moment,
        angles = copy(floe.angles),
        x_subfloe_points = copy(floe.x_subfloe_points),
        y_subfloe_points = copy(floe.y_subfloe_points),
        α = floe.α,
        u = floe.u,
        v = floe.v,
        ξ = floe.ξ,
        status = Status(floe.status.tag, copy(floe.status.fuse_idx)),
        id = floe.id,
        ghost_id = floe.ghost_id,
        parent_ids = copy(floe.parent_ids),
        ghosts = copy(floe.ghosts),
        fxOA= floe.fxOA,
        fyOA = floe.fyOA,
        trqOA = floe.trqOA,
        hflx_factor = floe.hflx_factor,
        overarea = floe.overarea,
        collision_force = copy(floe.collision_force),
        collision_trq = floe.collision_trq,
        stress_accum = copy(floe.stress_accum),
        stress_instant = copy(floe.stress_instant),
        strain = copy(floe.strain),
        p_dxdt = floe.p_dxdt,
        p_dydt = floe.p_dydt,
        p_dudt = floe.p_dudt,
        p_dvdt = floe.p_dvdt,
        p_dξdt = floe.p_dξdt,
        p_dαdt = floe.p_dαdt,
    )
    return f
end

"""
    translate!(coords, Δx, Δy)

Make a copy of given coordinates and translate by given deltas. 
Inputs:
    coords PolyVec{Float}
    vec <Vector{Real}>
Output:
    Updates given coords
"""
function translate(coords::PolyVec{FT}, Δx, Δy) where {FT<:AbstractFloat}
    new_coords = [[Vector{Float64}(undef, 2) for _ in eachindex(coords[1])]]
    for i in eachindex(coords[1])
        new_coords[1][i][1] = coords[1][i][1] + Δx
        new_coords[1][i][2] = coords[1][i][2] + Δy 
    end
    return new_coords
end

"""
    translate!(coords, Δx, Δy)

Translate each of the given coodinates by given deltas in place
Inputs:
    coords PolyVec{Float}
    vec <Vector{Real}>
Output:
    Updates given coords
"""
function translate!(coords::PolyVec{FT}, Δx, Δy) where {FT<:AbstractFloat}
    for i in eachindex(coords)
        for j in eachindex(coords[i])
            coords[i][j][1] += Δx
            coords[i][j][2] += Δy
        end
    end
    return
end

"""
    rotate_radians!(coords::PolyVec, α)

Rotate a polygon's coordinates by α radians around the origin.
Inputs:
    coords  <PolyVec{AbstractFloat}> polygon coordinates
    α       <Real> radians to rotate the coordinates
Outputs:
    Updates coordinates in place
"""
function rotate_radians!(coords::PolyVec, α)
    for i in eachindex(coords)
        for j in eachindex(coords[i])
            x, y = coords[i][j]
            coords[i][j][1] = cos(α)*x - sin(α)*y
            coords[i][j][2] = sin(α)*x + cos(α)*y
        end
    end
    return
end

"""
    hashole(coords::PolyVec{FT})

Determine if polygon coordinates have one or more holes
Inputs:
    coords <PolyVec{Float}>
Outputs:
    <Bool>
"""
function hashole(coords::PolyVec{FT}) where FT<:AbstractFloat
    return length(coords) > 1
end

"""
    hashole(poly::Polys)

Determine if polygon has one or more holes
Inputs:
    poly <Polygon> polygon
Outputs:
    <Bool> true if there is a hole in the polygons, else false
"""
function hashole(poly::Polys)
    return GI.nhole(poly) > 0
end 


function rmholes!(coords::PolyVec{FT}) where {FT<:AbstractFloat}
    if length(coords) > 1
        deleteat!(coords, 2:length(coords))
    end
end

function rmholes!(poly::Polys)
    deleteat!(poly.geom, 2:GI.nring(poly))
    return
end

#=
    _calc_moment_inertia(::Type{T} poly, cent, h; ρi = 920.0)

Calculate the mass moment of intertia from a polygon given the polygon, its centroid,
height, and the density of ice in the simulation. Answer will be of given type T.

Note: Assumes that first and last point within the coordinates are the same and will not
produce correct answer otherwise.

Based on paper: Marin, Joaquin."Computing columns, footings and gates through
moments of area." Computers & Structures 18.2 (1984): 343-349.
=#
function _calc_moment_inertia(
    ::Type{T},
    poly,
    cent,
    height;
    ρi = 920.0,
) where T
    xc, yc = GO._tuple_point(cent, T)
    Ixx, Iyy = zero(T), zero(T)
    x1, y1 = zero(T), zero(T)
    for (i, p2) in enumerate(GI.getpoint(GI.getexterior(poly)))
        (x2, y2) = GO._tuple_point(p2, T)
        x2, y2 = x2 - xc, y2 - yc
        if i == 1
            x1, y1 = x2, y2 
            continue
        end
        wi = (x1 - xc) * (y2 - yc) - (x2 - xc) * (y1 - yc)
        Ixx += wi * (y1^2 + y1 * y2 + y2^2)
        Iyy += wi * (x1^2 + x1 * x2 + x2^2)
        x1, y1 = x2, y2 
    end
    Ixx *= 1/12
    Iyy *= 1/12
    return abs(Ixx + Iyy) * T(height) * T(ρi) 
end

# Find the length of the maximum radius of a given polygon
function calc_max_radius(poly, cent, ::Type{T}) where T
    max_rad_sqrd = zero(T)
    Δx, Δy = GO._tuple_point(cent, T)
    for pt in GI.getpoint(GI.getexterior(poly))
        x, y = GO._tuple_point(pt, T)
        x, y = x - Δx, y - Δy
        rad_sqrd = x^2 + y^2
        if rad_sqrd > max_rad_sqrd
            max_rad_sqrd = rad_sqrd
        end
    end
    return sqrt(max_rad_sqrd)
end

"""
    which_vertices_match_points(ipoints, coords, atol)

Find which vertices in coords match given points
Inputs:
    points <Vector{Tuple{Float, Float} or Vector{Vector{Float}}}> points to
                match to vertices within polygon
    region  <Polygon> polygon 
    atol    <Float> distance vertex can be away from target point before being
                classified as different points
Output:
    Vector{Int} indices of points in polygon that match the intersection points
Note: 
    If last coordinate is a repeat of first coordinate, last coordinate index is
    NOT recorded.
"""
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

"""
euclidian_dist(c, idx2, idx1)

Calculate euclidean distance between two points within given coordinates
"""
euclidian_dist(c, idx2, idx1) = sqrt(
    (c[1][idx2][1] - c[1][idx1][1])^2 +
    (c[1][idx2][2] - c[1][idx1][2])^2 
)

"""
    which_points_on_edges(points, coords; atol = 1e-1)

Find which points are on the coordinates of the given polygon.
Inputs:
    points <Vector{Tuple{Float, Float} or Vector{Vector{Float}}}> points to
        match to edges within polygon
    coords  <PolVec> polygon coordinates
    atol    <Float> distance target point can be from an edge before being
                classified as not on the edge
"""
function which_points_on_edges(points, coords; atol = 1e-1)
    idxs = Vector{Int}()
    nedges = length(coords[1]) - 1
    npoints = length(points)
    if points[1] == points[end]
        npoints -= 1
    end
    for i in 1:nedges
        x1, y1 = coords[1][i]
        x2, y2 = coords[1][i+1]
        Δx_edge = x2 - x1
        Δy_edge = y2 - y1
        for j in 1:npoints
            xp, yp = points[j]
            Δx_point = xp - x1
            Δy_point = yp - y1
            if (!(j in idxs) && (
                (  # vertical edge
                    isapprox(Δx_edge, 0, atol = atol) &&
                    isapprox(Δx_point, 0, atol = atol) &&
                    0 < Δy_point / Δy_edge < 1
                ) ||
                (  # horizontal edge
                    isapprox(Δy_edge, 0, atol = atol) &&
                    isapprox(Δy_point, 0, atol = atol)
                    && 0 < Δx_point / Δx_edge < 1) ||
                (  # point is a vertex
                    isapprox(Δx_point, 0, atol = atol) &&
                    isapprox(Δy_point, 0, atol = atol)
                ) ||
                ( # has the same slope and is between edge points
                    isapprox(Δy_edge/Δx_edge, Δy_point/Δx_point, atol = atol) &&
                    0 < Δx_point / Δx_edge < 1 && 0 < Δy_point / Δy_edge < 1
                )
            ))
                push!(idxs, j)
            end
        end
    end
    sort!(idxs)
    return idxs
end

"""
    check_for_edge_mid(c, start, stop, shared_idx, shared_dist, running_dist)

Check if indices from start to stop index of given coords includes midpoint
given the shared distance and return midpoint if it exists in given range.
Inputs:
    c               <PolyVec> floe coordinates
    start           <Int> index of shared_index list to start search from
    stop            <Int> index of shared_index list to stop search at
    shared_idx      <Vector{Int}> list of indices of c used to calculate midpoint
    shared_dist     <Float> total length of edges considered from shared_idx
    running_dist    <Float> total length of edges traveled along in midpoint
                        search so far
Outputs:
    mid_x           <Float> x-coordinate of midpoint, Inf if midpoint not in
                        given range
    mid_y           <Float> y-coordinate of midpoint, Inf if midpoint not in
                        given range
    running_dist    <Float> sum of distances travelled along shared edges
"""
function check_for_edge_mid(c, start, stop, shared_idx, shared_dist,
    running_dist::FT,
) where FT
    mid_x = FT(Inf)
    mid_y = FT(Inf)
    for i in start:(stop-1)
        # Indices of c that are endpoints of current edge
        idx1 = shared_idx[i]
        idx2 = shared_idx[i + 1]
        # Lenght of edge
        edge_dist = euclidian_dist(c, idx2, idx1)
        # if midpoint is on current edge
        if running_dist + edge_dist >= shared_dist / 2
            frac = ((shared_dist / 2) - running_dist) / edge_dist
            mid_x = c[1][idx1][1] + (c[1][idx2][1] - c[1][idx1][1]) * frac
            mid_y = c[1][idx1][2] + (c[1][idx2][2] - c[1][idx1][2]) * frac
            break
        else  # move on to the next edge
            running_dist += edge_dist
        end
    end
    return mid_x, mid_y, running_dist
end

"""
    find_shared_edges_midpoint(c1, c2)

Find "midpoint" of shared polygon edges by distance
Inputs:
    c1      <PolVec> polygon coordinates for floe 1
    c2      <PolVec> polygon coordinates for floe 2
Outputs:
    mid_x   <Float> x-coordinate of midpoint
    mid_y   <Float> y-coordinate of midpoint
"""
function find_shared_edges_midpoint(c1::PolyVec{FT}, c2; atol = 1e-1) where {FT}
    # Find which points of c1 are on edges of c2
    shared_idx = which_points_on_edges(
        c1[1],
        c2;
        atol,
    )
    if shared_idx[1] == 1
         # due to repeated first point/last point
        push!(shared_idx, length(c1[1]))
    end
    shared_dist = FT(0)
    gap_idx = 1
    # Determine total length of edges that are shared between the two floes
    nshared_points = length(shared_idx)
    for i in 1:(nshared_points - 1)
        idx1 = shared_idx[i]
        idx2 = shared_idx[i + 1]
        if idx2 - idx1 == 1
            shared_dist += euclidian_dist(c1, idx2, idx1)
        elseif shared_dist > 0
            gap_idx = i + 1
            # Add distance wrapping around from last index to first
            shared_dist += euclidian_dist(c1, shared_idx[end], shared_idx[1])
        end
    end
    # Determine mid-point of shared edges by distance
    running_dist = FT(0)
    mid_x, mid_y, running_dist = check_for_edge_mid(c1, gap_idx, nshared_points,
        shared_idx, shared_dist, running_dist,
    )
    # Note that this assumes first and last point are the same 
    if isinf(mid_x)
        mid_x, mid_y, running_dist = check_for_edge_mid(c1, 1, gap_idx - 1,
            shared_idx, shared_dist, running_dist,
        )
    end
    return mid_x, mid_y
end