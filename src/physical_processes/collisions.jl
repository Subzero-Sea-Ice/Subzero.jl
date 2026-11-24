# Functions needed for collisions between floes, boundaries, and topography

"""
    calc_normal_force(...)

Calculate normal force for collision between polygons p1 and p2 given an overlapping region,
the area of that region, their intersection points in the region, and a force factor. 

## _Positional arguments_
- `p1::Polys`: first polygon
- `p2::Polys`: second polygon
- `region::Polys`: polygon for one region of intersection between p1 and p2
- `area::Float`: area of region
- `ipoints::Vector{Tuple{Float, Float}`: Points of intersection between p1 and p2
- `force_factor::Float`: Spring constant equivalent for collisions

##  _Returns_
- `Float`: normal force of collision
- `Float`: mean length of distance between intersection points
"""
function calc_normal_force(p1, p2, region, area, ipoints, force_factor::FT) where {FT<:AbstractFloat}
    force_dir = zeros(FT, 2)
    # Identify which region coordinates are the intersection points (ipoints)
    p = which_vertices_match_points(ipoints, region)
    m = length(p)
    # Calculate force direction
    Δl = FT(0)
    if m == 2  # Two intersection points
        idx1, idx2 = p
        pt1 = GI.getpoint(GI.getexterior(region), idx1)
        pt2 = GI.getpoint(GI.getexterior(region), idx2)
        Δx = FT(GI.x(pt2) - GI.x(pt1))
        Δy = FT(GI.y(pt2) - GI.y(pt1))
        Δl = sqrt(Δx^2 + Δy^2)
        if Δl > 0.1  # If overlap is large enough to consider
            force_dir .= [-Δy/Δl; Δx/Δl]
        end
    elseif m != 0  # Unusual number of intersection points
        Δl = _many_intersect_normal_force!(FT, force_dir, region, p1, force_factor)
    end
    # Check if direction of the force desceases overlap, else negate direction
    if Δl > 0.1
        p1new = _translate_poly(FT, p1, force_dir[1], force_dir[2])::Polys{FT}
        # Floe/boudary intersection after being moved in force direction
        new_regions_list = intersect_polys(p1new, p2, FT)
        # See if the area of overlap has increased in corresponding region
        for new_region in new_regions_list
            if check_intersects(new_region, region) && area_poly(new_region, FT)/area > 1
                force_dir .*= -1 
            end
        end
    end
    return force_dir * area * force_factor, Δl
end

#=
    _many_intersect_normal_force!(::Type{T}, force_dir, region, poly, force_factor)
    
Calculate the force direction (`force_dir`) given more than two points of intersection
within the region between two polygons (the first of which is `poly`).
=#
function _many_intersect_normal_force!(::Type{T}, force_dir, region, poly, force_factor) where T
    x1, y1 = zero(T), zero(T)
    Δl, n_pts = zero(T), 0
    Fn_tot = (zero(T), zero(T))
    # Loop over each edge within the overlap region
    for (i, p) in enumerate(GI.getpoint(GI.getexterior(region)))
        x2, y2 = T(GI.x(p)), T(GI.y(p))
        if i == 1
            x1, y1 = x2, y2
            continue
        end
        # Find the edge midpoint and calculate distance to first polygon
        xmid, ymid = 0.5 * (x2 + x1), 0.5 * (y2 + y1)
        dist = abs(dist_to_poly((xmid, ymid), poly, T))
        if dist < 1e-8  # Only consider region edge points on first polygon
            Δx, Δy = x2 - x1, y2 - y1
            mag = sqrt(Δx^2 + Δy^2)
            #= If force would push edge points closer to second polygon (past and out of the
            overlap region), switch the force direction =#
            xt = xmid + (-Δy / 100mag)
            yt = ymid + (Δx / 100mag)
            in_region = coveredby_poly((xt, yt), region) 
            f_sign = in_region ? 1 : -1
            # Calculate force from given edge and incorporate it into the total forces
            Fn = (f_sign * force_factor) .* (-Δy, Δx)
            Δl += mag
            n_pts += 1
            Fn_tot = Fn_tot .+ Fn
        end
        x1, y1 = x2, y2
    end
    # Take the average of the summed values
    if 0 < n_pts < (GI.npoint(region) - 1)
        Δl /= n_pts
        if Δl > 0.1  # If overlap is large enough to consider, set new force direction
            norm_Fn = sqrt(GI.x(Fn_tot)^2 + GI.y(Fn_tot)^2)
            force_dir[1] = GI.x(Fn_tot) / norm_Fn
            force_dir[2] = GI.y(Fn_tot) / norm_Fn
        end
    end
    return Δl
end

"""
    calc_elastic_forces(...)

Calculate normal forces, the point the force is applied, and the overlap area of
regions created from floe collisions 

## _Positional arguments_
- `p1::Polys`: first polygon
- `p2::Polys`: second polygon
- `regions::Vector{Polys}`: list of polygons for regions of intersection between p1 and p2
- `region_areas::Vector{Float}`: areas of regions
- `force_factor::Float`: Spring constant equivalent for collisions

##  _Returns_
- `force::Array{Float, n, 2}`: normal forces on each of the n regions greater than a minimum area
- `fpoint::Array{Float, n, 2}`: point force is applied on each of the n regions greater than a minimum area
- `overlap::Array{Float, n, 1}`: area of each of the n regions greater than a minimum area
- `Δl::Float`: mean length of distance between intersection points
"""
function calc_elastic_forces(p1, p2, regions, region_areas, force_factor::FT) where {FT<:AbstractFloat}
    ipoints = get_intersection_points(p1, p2)
    ncontact = 0
    if !isempty(ipoints) && length(ipoints) >= 2
        # Find overlapping regions greater than minumum area
        n1, n2 = GI.npoint(p1) - 1, GI.npoint(p2) - 1
        min_area = min(n1, n2) * 100 / 1.75
        for i in reverse(eachindex(region_areas))
            if region_areas[i] < min_area
                deleteat!(region_areas, i)
                deleteat!(regions, i)
            else
                ncontact += 1
            end
        end
    end
    # Calculate forces for each remaining region
    force = zeros(FT, ncontact, 2)
    fpoint = zeros(FT, ncontact, 2)
    overlap = ncontact > 0 ? region_areas : zeros(FT, ncontact)
    Δl_lst = zeros(FT, ncontact)
    for k in 1:ncontact
        if region_areas[k] != 0
            cx, cy = centroid_poly(regions[k], FT)
            fpoint[k, 1] = cx
            fpoint[k, 2] = cy
            normal_force, Δl = calc_normal_force(p1, p2, regions[k], region_areas[k], ipoints, force_factor)
            force[k, 1] = normal_force[1]
            force[k, 2] = normal_force[2]
            Δl_lst[k] = Δl
        end
    end
    return force, fpoint, overlap, Δl_lst
end

# Get velocity of a point, assumed to be on given floe. Includes the rotational velocity
# broken down into x- and y-directions.
function _get_velocity(
    floe::FloeType{FT},
    x::FT,
    y::FT,
) where {FT}
    u = floe.u + floe.ξ * (x - floe.centroid[1])
    v = floe.v + floe.ξ * (y - floe.centroid[2])
    return u, v
end

"""
    calc_friction_forces(...)

Calculate frictional force for collision between two floes or a floe and a
domain element.

## _Positional arguments_
- `ifloe::Floe`: first floe in collsion
- `jfloe::Union{Floe, DomainElement}`: either second floe or topography element/boundary element
- `fpoints::Array{Float, N, 2}`: x,y-coordinates of the point the force is applied on floe overlap region
- `normal::Array{Float, N, 2}`: x,y normal force applied on fpoint on floe overlap region
- `Δl::Vector`: mean length of distance between intersection points
- $CONSTS_DEF
- $ΔT_DEF

##  _Returns_
- `force::Array{Float, N, 2}`: frictional/tangential force of the collision in 
    x and y (each row) for each collision (each column)
"""
function calc_friction_forces(
    ifloe,
    jfloe,  # could be a boundary or a topography
    fpoints,
    normal::Matrix{FT},
    Δl,
    consts,
    Δt,
) where {FT}
    force = zeros(FT, size(fpoints, 1), 2)
    G = consts.E/(2*(1+consts.ν))
    for i in axes(fpoints, 1)
        px = fpoints[i, 1]
        py = fpoints[i, 2]
        nnorm = sqrt(normal[i, 1]^2 + normal[i, 2]^2)
        iu, iv = _get_velocity(ifloe, px, py)
        ju, jv = _get_velocity(jfloe, px, py)

        udiff = iu - ju
        vdiff = iv - jv

        vnorm = sqrt(udiff^2 + vdiff^2)
        xdir = FT(0)
        ydir = FT(0)
        if udiff != 0 || vdiff != 0
            xdir = udiff/vnorm
            ydir = vdiff/vnorm
        end
        dot_dir = xdir * udiff + ydir * vdiff
        xfriction = G *  Δl[i] * Δt * nnorm * xdir * -dot_dir
        yfriction = G *  Δl[i] * Δt * nnorm * ydir * -dot_dir
        norm_fric = sqrt(xfriction^2 + yfriction^2)
        if norm_fric > consts.μ * nnorm
            xfriction = -consts.μ * nnorm * xdir
            yfriction = -consts.μ * nnorm * ydir
        end
        force[i, 1] = xfriction
        force[i, 2] = yfriction
    end
    return force
end

"""
    add_interactions!(...)

Add interactions in-place to Floe's `interactions` field based on the number of interactions (`np`)
it had with other floe's/domain elements.

The interactions field has the following form for each interaction:

    `[Inf, xforce, yforce, xfpoints, yfpoints, overlaps]`

where the xforce and yforce are the forces, xfpoints and yfpoints are the location of the force
and overlaps is the overlap between the floe and boundary. The overlaps field is also added to
the floe's overarea field that describes the total overlapping area at any timestep. 

## _Positional arguments_
- `np::Int`: number of points to add to the interactions field
- `ifloe::Floe`: floe to add interactions to
- `idx::Int`: 
- `forces::Array{Float, n, 2}`: x and y-forces for each interaction
- `points::Array{Float, n, 2}`: x and y-location of interaction point
- `overs::Array{Float, n, 1}`:overlapping area between floes for each interaction
"""
function add_interactions!(np, ifloe, idx::FT, forces, points, overs) where FT
    inter_spots = size(ifloe.interactions, 1) - ifloe.num_inters
    @views for i in 1:np
        if forces[i, 1] != 0 || forces[i, 2] != 0
            ifloe.num_inters += 1
            if inter_spots < 1
                ifloe.interactions = vcat(
                    ifloe.interactions,
                    zeros(FT, np - i + 1 , 7)
                )
                inter_spots += np - i
            end
            ifloe.interactions[ifloe.num_inters, floeidx] = idx
            ifloe.interactions[ifloe.num_inters, xforce:yforce] .=
                forces[i, :]
            ifloe.interactions[ifloe.num_inters, xpoint:ypoint] .=
                points[i, :]
            ifloe.interactions[ifloe.num_inters, torque] = FT(0)
            ifloe.interactions[ifloe.num_inters, overlap] = overs[i]
            ifloe.overarea += overs[i]
            inter_spots -= 1
        end
    end
    return
end

"""
    floe_floe_interaction!(...)

If the two floes interact, update floe i's interactions accordingly. Floe j is
not update here so that the function can be parallelized.

## _Positional arguments_
- `ifloe::Floe`: first floe in potential interaction
- `i::Int`: index of ifloe in model's list of floes 
- `jfloe::Floe`: second floe in potential interaction
- `j::Int`: index of jfloe in model's list of floes 
- `nfloes::Int`: number of non-ghost floes in the simulation this timestep
- $CONSTS_DEF
- $ΔT_DEF
- `max_overlap::Float`: Percent two floes can overlap before marking them for fusion
"""
function floe_floe_interaction!(
    ifloe,
    i,
    jfloe,
    j,
    consts,
    Δt,
    max_overlap::FT,
) where {FT<:AbstractFloat}
    inter_regions = intersect_polys(ifloe.poly, jfloe.poly, FT)
    region_areas = Vector{FT}(undef, length(inter_regions))
    total_area = FT(0)
    for i in eachindex(inter_regions)
        a = FT(area_poly(inter_regions[i], FT))
        region_areas[i] = a
        total_area += a
    end
    if total_area > 0
        # Floes overlap too much - mark floes to be fused
        if max(total_area/ifloe.area, total_area/jfloe.area) > max_overlap
            ifloe.status.tag = fuse
            push!(ifloe.status.fuse_idx, j)
        else
            # Constant needed to calculate force
            ih = ifloe.height
            ir = FT(sqrt(ifloe.area))
            jh = jfloe.height
            jr = FT(sqrt(jfloe.area))
            force_factor = if ir>1e5 || jr>1e5
                consts.E*min(ih, jh)/min(ir, jr)
            else
                consts.E*(ih*jh)/(ih*jr+jh*ir)
            end
            # Calculate normal forces, force points, and overlap areas
            normal_forces, fpoints, overlaps, Δl = calc_elastic_forces(
                ifloe.poly,
                jfloe.poly,
                inter_regions,
                region_areas,
                force_factor
            )
            #= Calculate frictional forces at each force point - based on
            velocities at force points =#
            np = size(fpoints, 1)
            if np > 0
                friction_forces = calc_friction_forces(
                    ifloe,
                    jfloe,
                    fpoints,
                    normal_forces,
                    Δl,
                    consts,
                    Δt,
                )
                # Calculate total forces and update ifloe's interactions
                forces = normal_forces .+ friction_forces
                add_interactions!(np, ifloe, j, forces, fpoints, overlaps)
            end
        end
    end
    return
end

#=
Collision forces between pairs of floes in a `FixedWidthFloes`, one pair at a time, so that a
GPU kernel can compute many pairs at once. The functions write their results to buffers with
one row per pair (`PairForces`), using scratch space with one row per pair (a
`CollisionWorkspace`), instead of adding them to the floes' interactions. They take the
polygons as `GOG.ArrayPolygon`s, so that they work for floes and domain elements alike.
=#

"""
    CollisionWorkspace(FT, nrows, max_points; max_ipoints = 16, max_regions = 4)
    CollisionWorkspace(workspace, k)

Scratch space for the collision forces between pairs of polygons of float type `FT` with up
to `max_points` points each (for floes, `size(floes.poly, 2)`), with `nrows` rows: one per
pair that is processed at the same time. The first form makes the arrays on the CPU (move
them to a device with `adapt`); the second refers to row `k` of the same arrays, e.g. the
row of one pair in a kernel.

Each row has room for `max_ipoints` points where the polygons' edges intersect and
`max_regions` overlap regions. The regions have room for `2max_points + max_ipoints + 1`
points each and the node lists for clipping for `max_points + max_ipoints` nodes, which is
enough as long as the intersection points fit. The functions that use the workspace return
`GOG.CLIP_TOO_SMALL` if a pair needs more room.

Fields:
- `lists`: node lists for clipping, used by both clippings in turn
- `regions`: the overlap regions of the two polygons
- `new_regions`: the overlap regions after the first polygon is moved in the direction of a
  normal force, in `calc_normal_force`
- `ipoints`: the points where the polygons' edges intersect (`ipoints[k, :, :]`)
- `vertex_idx`: the vertices of a region that are such points (`vertex_idx[k, :]`)
- `moved`: the first polygon, moved in the direction of a normal force (`moved[k, :, :]`)
- `Δl`: for each force, the mean length between the intersection points of its region
  (`Δl[k, r]`), which the friction needs
"""
struct CollisionWorkspace{L<:GOG.ClipLists, R<:GOG.PolyRegions, A3<:AbstractArray{<:Any, 3},
        I2<:AbstractMatrix{Int32}, A2<:AbstractMatrix}
    lists::L
    regions::R
    new_regions::R
    ipoints::A3
    vertex_idx::I2
    moved::A3
    Δl::A2
    k::Int
end

Adapt.@adapt_structure CollisionWorkspace

function CollisionWorkspace(::Type{FT}, nrows, max_points; max_ipoints = 16,
        max_regions = 4) where FT
    region_capacity = 2max_points + max_ipoints + 1
    return CollisionWorkspace(
        GOG.ClipLists(FT, nrows, max_points + max_ipoints),
        GOG.PolyRegions(FT, nrows, max_regions, region_capacity),
        GOG.PolyRegions(FT, nrows, max_regions, region_capacity),
        zeros(FT, nrows, max_ipoints, 2),
        zeros(Int32, nrows, max_ipoints),
        zeros(FT, nrows, max_points, 2),
        zeros(FT, nrows, max_regions),
        0,
    )
end

CollisionWorkspace(w::CollisionWorkspace, k) = CollisionWorkspace(GOG.ClipLists(w.lists, k),
    GOG.PolyRegions(w.regions, k), GOG.PolyRegions(w.new_regions, k), w.ipoints, w.vertex_idx,
    w.moved, w.Δl, Int(k))

"""
    PairForces(FT, nrows, max_forces)
    PairForces(forces, k)

Collision forces between pairs of polygons, with one row per pair and room for `max_forces`
forces per pair: force `r` of the pair in row `k` is `force[k, r, :]`, it acts at the point
`point[k, r, :]`, and it comes from an overlap region with area `overlap[k, r]`. Forces that
are zero are included. The first form makes the arrays on the CPU (move them to a device
with `adapt`); the second refers to row `k` of the same arrays.
"""
struct PairForces{A3<:AbstractArray{<:Any, 3}, A2<:AbstractMatrix}
    force::A3
    point::A3
    overlap::A2
    k::Int
end

Adapt.@adapt_structure PairForces

PairForces(::Type{FT}, nrows, max_forces) where FT = PairForces(
    zeros(FT, nrows, max_forces, 2), zeros(FT, nrows, max_forces, 2),
    zeros(FT, nrows, max_forces), 0)

PairForces(f::PairForces, k) = PairForces(f.force, f.point, f.overlap, Int(k))

# The number of forces per pair that `forces` has room for
max_forces(f::PairForces) = size(f.force, 2)

#=
Velocity of the point (x, y) of floe i, including the floe's rotation.

TODO: the rotation's part looks wrong. A floe that rotates with angular velocity ξ around its
centroid (cx, cy) moves the point with (-ξ(y - cy), ξ(x - cx)), not (ξ(x - cx), ξ(y - cy)).
=#
function _get_velocity(floes::FixedWidthFloes, i, x, y)
    u = floes.u[i] + floes.ξ[i] * (x - floes.centroid[i, 1])
    v = floes.v[i] + floes.ξ[i] * (y - floes.centroid[i, 2])
    return u, v
end

#=
Friction between two bodies (floes, or a floe and a domain element) that press against each
other with the normal force (nx, ny), at a point where they move with velocities (iu, iv) and
(ju, jv), when the mean length between the intersection points of their overlap region is Δl.
It opposes the relative velocity, and is at most the coefficient of friction times the normal
force.
=#
function _calc_friction_force(iu, iv, ju, jv, nx::FT, ny::FT, Δl, consts, Δt) where FT
    G = consts.E/(2*(1+consts.ν))
    nnorm = sqrt(nx^2 + ny^2)
    udiff = iu - ju
    vdiff = iv - jv
    vnorm = sqrt(udiff^2 + vdiff^2)
    xdir = FT(0)
    ydir = FT(0)
    if udiff != 0 || vdiff != 0
        xdir = udiff/vnorm
        ydir = vdiff/vnorm
    end
    dot_dir = xdir * udiff + ydir * vdiff
    xfriction = G * Δl * Δt * nnorm * xdir * -dot_dir
    yfriction = G * Δl * Δt * nnorm * ydir * -dot_dir
    norm_fric = sqrt(xfriction^2 + yfriction^2)
    if norm_fric > consts.μ * nnorm
        xfriction = -consts.μ * nnorm * xdir
        yfriction = -consts.μ * nnorm * ydir
    end
    return xfriction, yfriction
end

#=
The direction of the normal force on the overlap region `region` of `poly` and another
polygon, for a region with other than two intersection points: the sum of the normals of the
region's edges that are on `poly`, each pointing into the region and as long as its edge,
normalised. Returns the mean length of those edges, Δl, and the direction, which is
`force_dir` unless Δl > 0.1.

TODO: this has three problems, to be fixed:
- The edges on `poly` are those whose midpoints are less than 1e-8 from it. That tolerance is
  absolute, but the rounding errors of the distance grow with the coordinates: with Float32
  floes no edge passes, so the force is zero, and with Float64 floes edges start to be missed
  from coordinates of about 3e7 m.
- Which way an edge's normal points is decided at a point 1 cm from the edge, which is lost
  to rounding with Float32 coordinates from about 1e5 m.
- If the normals cancel out, e.g. for a narrow floe across a wider one, the direction is NaN
  (0 / 0). `calc_normal_force` then returns `GOG.CLIP_FAILED`.
=#
function _many_intersect_normal_force(::Type{T}, force_dir,
        region::GOG.ArrayPolygon{<:AbstractArray{T}}, poly::GOG.ArrayPolygon{<:AbstractArray{T}},
        force_factor) where T
    x1, y1 = zero(T), zero(T)
    Δl, n_pts = zero(T), 0
    Fn_tot = (zero(T), zero(T))
    # Loop over each edge of the overlap region
    for i in 1:GOG.npoint(region)
        x2, y2 = GOG.getpoint(region, i)
        if i == 1
            x1, y1 = x2, y2
            continue
        end
        # The edge's midpoint and its distance to `poly`
        xmid, ymid = T(0.5) * (x2 + x1), T(0.5) * (y2 + y1)
        dist = abs(dist_to_poly((xmid, ymid), poly, T))
        if dist < T(1e-8)  # only the region's edges on `poly`
            Δx, Δy = x2 - x1, y2 - y1
            mag = sqrt(Δx^2 + Δy^2)
            # Reverse the normal if it points out of the region
            xt = xmid + (-Δy / 100mag)
            yt = ymid + (Δx / 100mag)
            in_region = coveredby_poly((xt, yt), region)
            f_sign = in_region ? 1 : -1
            s = f_sign * force_factor
            Fn = (s * -Δy, s * Δx)
            Δl += mag
            n_pts += 1
            Fn_tot = (Fn_tot[1] + Fn[1], Fn_tot[2] + Fn[2])
        end
        x1, y1 = x2, y2
    end
    # The mean of the summed values
    if 0 < n_pts < (GOG.npoint(region) - 1)
        Δl /= n_pts
        if Δl > T(0.1)  # large enough to set the direction
            norm_Fn = sqrt(Fn_tot[1]^2 + Fn_tot[2]^2)
            force_dir = (Fn_tot[1] / norm_Fn, Fn_tot[2] / norm_Fn)
        end
    end
    return Δl, force_dir
end

#=
The normal force on the overlap region `region` of polygons `p1` and `p2`, whose area is
`area`, given the `n_ipoints` points where the polygons' edges intersect, in
`ws.ipoints[ws.k, :, :]`. It pushes p1 out of the region: its direction is perpendicular to
the line between the region's two intersection points, or comes from
`_many_intersect_normal_force` if the region has another number of them, and it's reversed if
moving p1 in that direction makes the overlap larger. Its size is `area * force_factor`. Uses
`ws` for the clipping.

Returns the force, the mean length between the intersection points Δl, and 0, or a negative
code as the third value: `GOG.CLIP_TOO_SMALL` if `ws` is too small, and `GOG.CLIP_FAILED` if
the clipping failed or the direction is NaN.

TODO: the direction is reversed once for each region of the overlap after the move that
intersects `region` and is larger, so two such regions cancel each other out.
=#
function calc_normal_force(ws::CollisionWorkspace, p1::GOG.ArrayPolygon,
        p2::GOG.ArrayPolygon, region::GOG.ArrayPolygon, area, n_ipoints,
        force_factor::FT) where FT
    k = ws.k
    force_dir = (zero(FT), zero(FT))
    # Which of the region's vertices are intersection points
    m = which_vertices_match_points!(ws.vertex_idx, k, ws.ipoints, n_ipoints, region)
    m < 0 && return force_dir, zero(FT), GOG.CLIP_TOO_SMALL
    Δl = zero(FT)
    if m == 2  # two intersection points
        x1, y1 = GOG.getpoint(region, ws.vertex_idx[k, 1])
        x2, y2 = GOG.getpoint(region, ws.vertex_idx[k, 2])
        Δx = x2 - x1
        Δy = y2 - y1
        Δl = sqrt(Δx^2 + Δy^2)
        if Δl > FT(0.1)  # if the overlap is large enough to consider
            force_dir = (-Δy/Δl, Δx/Δl)
        end
    elseif m != 0  # another number of intersection points
        Δl, force_dir = _many_intersect_normal_force(FT, force_dir, region, p1, force_factor)
    end
    # Reverse the direction if moving p1 in it makes the overlap larger
    if Δl > FT(0.1)
        # Moving p1 by NaN makes the clipping fail
        (isnan(force_dir[1]) || isnan(force_dir[2])) && return force_dir, Δl, GOG.CLIP_FAILED
        GOG.npoint(p1) > size(ws.moved, 2) && return force_dir, Δl, GOG.CLIP_TOO_SMALL
        p1new = _translate_poly!(GOG.ArrayPolygon(ws.moved, k, 0), p1, force_dir[1],
            force_dir[2])
        n_new = intersect_polys!(ws.new_regions, ws.lists, p1new, p2)
        n_new < 0 && return force_dir, Δl, n_new
        for r in 1:n_new
            new_region = GOG.region(ws.new_regions, r)
            if check_intersects(new_region, region) && area_poly(new_region, FT)/area > 1
                force_dir = (force_dir[1] * -1, force_dir[2] * -1)
            end
        end
    end
    force = (force_dir[1] * area * force_factor, force_dir[2] * area * force_factor)
    return force, Δl, 0
end

#=
The normal forces on the first `n_regions` overlap regions of polygons `p1` and `p2`, which
are in `ws.regions`. Each region that is large enough gets a force, which is written to row
`ws.k` of `out` with its point (the region's centroid) and the region's area; its Δl (see
`calc_normal_force`) is written to `ws.Δl`. Only if the polygons' edges intersect at two or
more points do the regions get forces.

Returns the number of forces, or a negative code (see `calc_normal_force`).
=#
function calc_elastic_forces!(out::PairForces, ws::CollisionWorkspace, p1::GOG.ArrayPolygon,
        p2::GOG.ArrayPolygon, n_regions, force_factor::FT) where FT
    k = ws.k
    n_ipoints = get_intersection_points!(ws.ipoints, k, p1, p2, FT)
    n_ipoints < 0 && return GOG.CLIP_TOO_SMALL
    n_ipoints >= 2 || return 0
    # Regions smaller than this get no force
    n1, n2 = GOG.npoint(p1) - 1, GOG.npoint(p2) - 1
    min_area = FT(min(n1, n2) * 100) / FT(1.75)
    n_forces = 0
    for r in 1:n_regions
        region = GOG.region(ws.regions, r)
        area = area_poly(region, FT)
        area < min_area && continue
        (n_forces == max_forces(out) || n_forces == size(ws.Δl, 2)) && return GOG.CLIP_TOO_SMALL
        n_forces += 1
        force, Δl = (zero(FT), zero(FT)), zero(FT)
        px, py = zero(FT), zero(FT)
        if area != 0
            px, py = centroid_poly(region, FT)
            force, Δl, err = calc_normal_force(ws, p1, p2, region, area, n_ipoints,
                force_factor)
            err < 0 && return err
        end
        out.force[k, n_forces, 1], out.force[k, n_forces, 2] = force
        out.point[k, n_forces, 1], out.point[k, n_forces, 2] = px, py
        out.overlap[k, n_forces] = area
        ws.Δl[k, n_forces] = Δl
    end
    return n_forces
end

"""
    floe_floe_interaction!(out::PairForces, ws::CollisionWorkspace, floes::FixedWidthFloes,
        i, j, consts, Δt, max_overlap)

The collision forces on floe i from floe j: a normal force and friction for each region
where they overlap enough. They are written to row `ws.k` of `out` (see `PairForces`),
using row `ws.k` of `ws` (see `CollisionWorkspace`), so that a kernel can compute many pairs
at once; floe i's interactions aren't changed. `max_overlap` is the fraction of either
floe's area that the floes can overlap before they should fuse.

Returns the number of forces n and whether the floes should fuse (then n is 0), or a
negative code and `false` if the workspace is too small (`GOG.CLIP_TOO_SMALL`) or the
clipping failed (`GOG.CLIP_FAILED`).
"""
function floe_floe_interaction!(out::PairForces, ws::CollisionWorkspace,
        floes::FixedWidthFloes{FT}, i, j, consts, Δt, max_overlap) where FT
    k = ws.k
    p1, p2 = _floe_poly(floes, i), _floe_poly(floes, j)
    n_regions = intersect_polys!(ws.regions, ws.lists, p1, p2)
    n_regions < 0 && return n_regions, false
    total_area = zero(FT)
    for r in 1:n_regions
        total_area += area_poly(GOG.region(ws.regions, r), FT)
    end
    total_area > 0 || return 0, false
    # The floes overlap too much: they should fuse
    if max(total_area/floes.area[i], total_area/floes.area[j]) > max_overlap
        return 0, true
    end
    # Constant needed to calculate the forces
    ih = floes.height[i]
    ir = sqrt(floes.area[i])
    jh = floes.height[j]
    jr = sqrt(floes.area[j])
    force_factor = if ir > FT(1e5) || jr > FT(1e5)
        consts.E*min(ih, jh)/min(ir, jr)
    else
        consts.E*(ih*jh)/(ih*jr+jh*ir)
    end
    n = calc_elastic_forces!(out, ws, p1, p2, n_regions, force_factor)
    n < 0 && return n, false
    # Add the friction at each force's point, from the floes' velocities there
    for r in 1:n
        px, py = out.point[k, r, 1], out.point[k, r, 2]
        nx, ny = out.force[k, r, 1], out.force[k, r, 2]
        iu, iv = _get_velocity(floes, i, px, py)
        ju, jv = _get_velocity(floes, j, px, py)
        fx, fy = _calc_friction_force(iu, iv, ju, jv, nx, ny, ws.Δl[k, r], consts, Δt)
        out.force[k, r, 1] = nx + fx
        out.force[k, r, 2] = ny + fy
    end
    return n, false
end

"""
    floe_domain_element_interaction!(floe, boundary::OpenBoundary, ...)

Dispatch off of `OpenBoundary`.

If given floe insersects with `OpenBoundary`, the floe is set to be removed
from the simulation.

## _Positional arguments_
- `floe::Floe`: floe interacting with boundary
- `boundary::OpenBoundary`: coordinates of boundary
- others inputs are NOT used in this dispatch
"""
function floe_domain_element_interaction!(
    floe::FloeType{FT},
    boundary::OpenBoundary,
    element_idx,
    consts,
    Δt,
    max_overlap,
) where FT
    # Check if the floe and boundary actually overlap
    FT_area_poly(p) = area_poly(p, FT)
    inter_area = sum(FT_area_poly, intersect_polys(floe.poly, boundary.poly, FT); init = 0.0)
    if inter_area > 0
        floe.status.tag = remove
    end
    return
end

"""
    floe_domain_element_interaction!(floe, boundary::PeriodicBoundary, ...)

Dispatch off of `PeriodicBoundary`.

If a given floe intersects with a periodic boundary, nothing happens at this
point. Periodic floes pass through boundaries using ghost floes.

## _Positional arguments_
- inputs are NOT used in this dispatch
"""
function floe_domain_element_interaction!(
    floe,
    ::PeriodicBoundary,
    element_idx,
    consts,
    Δt,
    max_overlap,
)
    return
end

"""
    floe_domain_element_interaction!(floe, element::Union{CollisionBoundary, MovingBoundary, TopographyElement}, ...)

Dispatch off of `CollisionBoundary`, `MovingBoundary`, or `TopographyElement`.

If floe intersects with given element (either collision boundary, moving boundary, or
topography element), floe interactions field is updated with details of the collision. 

## _Positional arguments_
- `floe::Floe`: floe interacting with element
- `element::Union{CollisionBoundary, MovingBoundary, TopographyElement}`: element floe is colliding with
- `elem_idx::Int`: negative integer value notating which type of element the floe is colliding with
    (-1 to -4 for NSEW walls and lower numbers for each topographic element)
- `consts::Constants`: model constants needed for calculations
- $ΔT_DEF
- `max_overlap::Float`: Percent a floe can overlap with an element before being killed/removed
"""
function floe_domain_element_interaction!(
    floe,
    element::Union{
        CollisionBoundary,
        MovingBoundary,
        TopographyElement,
    },
    elem_idx,
    consts,
    Δt,
    max_overlap::FT,
) where {FT}
    inter_regions = intersect_polys(floe.poly, element.poly, FT)
    region_areas = Vector{FT}(undef, length(inter_regions))
    max_area = FT(0)
    for i in eachindex(inter_regions)
        a = area_poly(inter_regions[i], FT)
        region_areas[i] = a
        if a > max_area
            max_area = a
        end
    end

    if max_area > 0
        # Regions overlap too much
        if maximum(region_areas)/floe.area > max_overlap
            floe.status.tag = remove
        else
            # Constant needed for force calculations
            force_factor = consts.E * floe.height / sqrt(floe.area)
            # Calculate normal forces, force points, and overlap areas
            normal_forces, fpoints, overlaps, Δl =  calc_elastic_forces(
                floe.poly,
                element.poly,
                inter_regions,
                region_areas,
                force_factor,
            )
            _normal_direction_correct!(normal_forces, fpoints, element)
            # Calculate frictional forces at each force point
            np = size(fpoints, 1)
            if np > 0
                friction_forces = calc_friction_forces(
                    floe,
                    element,
                    fpoints,
                    normal_forces,
                    Δl,
                    consts,
                    Δt,
                )
                # Calculate total forces and update ifloe's interactions
                forces = normal_forces .+ friction_forces
                add_interactions!(np, floe, elem_idx, forces, fpoints, overlaps)
            end
        end
    end
    return
end

"""
    update_boundaries!(domain)
    
Update each boundary in the domain. For now, this simply means updating the position of any `MovingBoundary` elements.
"""
function update_boundaries!(domain, Δt)
    _update_boundary!(domain.north, Δt)
    _update_boundary!(domain.south, Δt)
    _update_boundary!(domain.east, Δt)
    _update_boundary!(domain.west, Δt)
    return
end

"""
    floe_domain_interaction!(floe, domain::Domain, ...)

If the floe interacts with the domain, update the floe accordingly. Dispatches
on different boundary types and topography elements within the domain.

## _Positional arguments_
- `floe::Floe`: floe interacting with element
- $DOMAIN_DEF
- $CONSTS_DEF
- $ΔT_DEF
- `max_overlap::Float`: Percent a floe can overlap with an element before being killed/removed
"""
function floe_domain_interaction!(
    floe,
    domain::Domain,
    consts,
    Δt,
    max_overlap::FT,
) where {FT<:AbstractFloat}
    centroid = floe.centroid
    rmax = floe.rmax
    nbound = domain.north
    sbound = domain.south
    ebound = domain.east
    wbound = domain.west

    if centroid[2] + rmax > nbound.val
        floe_domain_element_interaction!(
            floe,
            nbound,
            -1,  # north wall index
            consts,
            Δt,
            max_overlap,
        )
    end
    if centroid[2] - rmax < sbound.val
        floe_domain_element_interaction!(
            floe,
            sbound,
            -2,  # south wall index
            consts,
            Δt,
            max_overlap,
        )
    end
    if centroid[1] + rmax > ebound.val
        floe_domain_element_interaction!(
            floe,
            ebound,
            -3,  # east wall index
            consts,
            Δt,
            max_overlap,
        )
    end
    if centroid[1] - rmax < wbound.val
        floe_domain_element_interaction!(
            floe,
            wbound,
            -4,  # west wall index
            consts,
            Δt,
            max_overlap,
        )
    end

    for (i, topo_element) in enumerate(domain.topography)
        if sum((topo_element.centroid .- floe.centroid).^2) < (topo_element.rmax + floe.rmax)^2
            floe_domain_element_interaction!(
                floe,
                topo_element,
                -(4 + i),  # topograph elements indices
                consts,
                Δt,
                max_overlap,
            )
        end
    end
    return
end

#Calculate a floe's torque based on the interactions.
function _calc_torque!(
    floe::Union{LazyRow{<:Floe{FT}}, Floe{FT}},
) where {FT<:AbstractFloat}
    inters = floe.interactions
    for i in 1:floe.num_inters
        xp = inters[i, xpoint] .- floe.centroid[1]
        yp = inters[i, ypoint] .- floe.centroid[2]
        xf = inters[i, xforce]
        yf = inters[i, yforce]
        # cross product of direction and force where third element of both is 0
        floe.interactions[i, torque] = xp * yf - yp * xf
    end
    return
end


#=
Determine if two floes could potentially interact using the two centroid and two
radii to form a bounding circle.
=#
_potential_interaction(centroid1, centroid2, rmax1, rmax2) =
    ((centroid1[1] - centroid2[1])^2 + (centroid1[2] - centroid2[2])^2) < (rmax1 + rmax2)^2

"""
    timestep_collisions!(...)

Resolves collisions between pairs of floes and domain elements and calculates the forces and torques
caused by those collisions.

## _Positional arguments_
- $FLOES_DEF
- `n_init_floes:Int`: number of floes without ghost floes
- $DOMAIN_DEF
- $CONSTS_DEF
- $ΔT_DEF
- `collision_settings::CollisionSettings`: simulation's collision settings
- `spinlock::Thread.SpinLock`: lock's multi-threading to prevent parallel errors
"""
function timestep_collisions!(
    floes::StructArray{<:Floe{FT}},
    n_init_floes,
    domain,
    consts,
    Δt,
    collision_settings,
    spinlock,
) where {FT<:AbstractFloat}
    collide_pairs = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()
    # floe-floe collisions for floes i and j where i<j
    Threads.@threads for i in eachindex(floes)
        # reset collision values
        fill!(floes.collision_force[i], FT(0))
        floes.collision_trq[i] = FT(0.0)
        floes.num_inters[i] = 0
        for j in i+1:length(floes)
            id_pair, ghost_id_pair =
                if floes.id[i] > floes.id[j]
                    (floes.id[i], floes.id[j]), (floes.ghost_id[i], floes.ghost_id[j])
                else
                    (floes.id[j], floes.id[i]), (floes.ghost_id[j], floes.ghost_id[i])
                end
            # If checking two distinct floes (i.e. not a parent ghost pair) and they are in close proximity continue
            if (id_pair[1] != id_pair[2]) && _potential_interaction(
                floes.centroid[i],
                floes.centroid[j],
                floes.rmax[i],
                floes.rmax[j],
            )
                # Add floes/ghosts to collision list if they didn't exist before
                (g1, g2) = Threads.lock(spinlock) do
                    get!(collide_pairs, id_pair, ghost_id_pair)
                end
                #=
                    New collision or floe and ghost colliding with same floe ->
                    not a repeat collision
                =#
                if (
                    (ghost_id_pair[1] == g1) && (ghost_id_pair[2] == g2) ||
                    (ghost_id_pair[1] == g1) ⊻ (ghost_id_pair[2] == g2)
                )
                    floe_floe_interaction!(
                        get_floe(floes, i),
                        i,
                        get_floe(floes, j),
                        j,
                        consts,
                        Δt,
                        collision_settings.floe_floe_max_overlap,
                    )
                end
            end
        end
        floe_domain_interaction!(
            get_floe(floes, i),
            domain,
            consts,
            Δt,
            collision_settings.floe_domain_max_overlap,
        )
    end
    # Move compression boundaries if they exist
    update_boundaries!(domain, Δt)
    # Update floes not directly calculated above where i>j - can't be parallel
    for i in eachindex(floes)
        # Update fuse information
        if floes.status[i].tag == fuse
            for idx in floes.status[i].fuse_idx
                floes.status[idx].tag = fuse
                push!(floes.status[idx].fuse_idx, i)
            end
        end
        # Update interaction information
        i_inters = floes.interactions[i]
        if floes.num_inters[i] > 0
            # Loop over each interaction with floe i
            for inter_idx in 1:floes.num_inters[i]
                j = i_inters[inter_idx, floeidx]  # Index of floe to update
                # not boundaries and hasn't been updated yet
                if j <= length(floes) && j > i
                    jidx = Int(j)
                    # add matching interaction (j with i)
                    add_interactions!(1, get_floe(floes, jidx), FT(i),
                        i_inters[inter_idx:inter_idx, xforce:yforce],
                        i_inters[inter_idx:inter_idx, xpoint:ypoint],
                        i_inters[inter_idx:inter_idx, overlap:overlap]
                    )
                    # equal and opposite forces
                    j_np = floes.num_inters[jidx]
                    floes.interactions[jidx][j_np, xforce:yforce] *= -1
                end
            end
        end
    end
    # Calculate total on parents by summing interactions on parent and children
    for i in 1:n_init_floes
        if !isempty(floes.ghosts[i])
            for g in floes.ghosts[i]
                gnp = floes.num_inters[g]
                # shift interaction points to parent locations
                floes.interactions[g][1:gnp, xpoint] .-=
                    floes.centroid[g][1] - floes.centroid[i][1]
                floes.interactions[g][1:gnp, ypoint] .-=
                    floes.centroid[g][2] - floes.centroid[i][2]
                # add interactions
                add_interactions!(gnp, get_floe(floes, i), i,
                    floes.interactions[g][1:gnp, xforce:yforce],
                    floes.interactions[g][1:gnp, xpoint:ypoint],
                    floes.interactions[g][1:gnp, overlap]
                )
                inp = floes.num_inters[i]
                # set interaction indices to correct values
                floes.interactions[i][(inp - gnp + 1):inp, floeidx] .= 
                    floes.interactions[g][1:gnp, floeidx]
            end
        end
        # calculate interactions torque and total forces / torque
        _calc_torque!(get_floe(floes, i))
        floes.collision_force[i][1] += sum(
            @view floes.interactions[i][1:floes.num_inters[i], xforce]
        )
        floes.collision_force[i][2] += sum(
            @view floes.interactions[i][1:floes.num_inters[i], yforce]
        )
        floes.collision_trq[i] += sum(
            @view floes.interactions[i][1:floes.num_inters[i], torque]
        )
    end
    return
end

#=
If the given element intersects with the boundary, add ghosts of the element and 
any of its existing ghosts to the `floes` list. The ghost will be translated from
the "parent" floe with `trans_vec`.
=#
function ghosts_on_bounds!(
    floes::FLT,
    elem_idx,
    boundary,
    trans_vec,
) where {FT <: AbstractFloat, FLT <: StructArray{<:Floe{FT}}}
    nfloes = length(floes)
    nghosts = 1
    if !isempty(intersect_polys(floes.poly[elem_idx], boundary.poly, FT))
        # ghosts of existing ghosts and original element
        for i in floes.ghosts[elem_idx]
            push!(floes, deepcopy_floe(get_floe(floes, i)))
            nghosts += 1
        end
        push!(floes, deepcopy_floe(get_floe(floes, elem_idx)))
        for i in (nfloes + 1):(nfloes + nghosts)
            _translate_floe!(FT, get_floe(floes, i), trans_vec...)
        end
    end
    return
end

#=
Find ghosts of given element and its known ghosts through an eastern or western
periodic boundary. If element's centroid isn't within the domain in the
north/south direction, swap it with its ghost since the ghost's centroid must
then be within the domain. `elem_idx`is the floe of floe of interest's index within
the floe list.
=#
function find_ghosts!(
    floes,
    elem_idx,
    ebound::PeriodicBoundary{East, FT},
    wbound::PeriodicBoundary{West, FT},
) where {FT <: AbstractFloat}
    Lx = ebound.val - wbound.val
    nfloes = length(floes)
    # passing through western boundary
    if (floes.centroid[elem_idx][1] - floes.rmax[elem_idx] < wbound.val)
        ghosts_on_bounds!(floes, elem_idx, wbound, [Lx, FT(0)])
    # passing through eastern boundary
    elseif (floes.centroid[elem_idx][1] + floes.rmax[elem_idx] > ebound.val)
        ghosts_on_bounds!(floes, elem_idx, ebound, [-Lx, FT(0)])
    end
    # if element centroid isn't in domain's east/west direction, swap with ghost
    new_nfloes = length(floes)
    if new_nfloes > nfloes
        if floes.centroid[elem_idx][1] < wbound.val
            _translate_floe!(FT, get_floe(floes, elem_idx), Lx, zero(FT))
            _translate_floe!(FT, get_floe(floes, new_nfloes), -Lx, zero(FT))
        elseif ebound.val < floes.centroid[elem_idx][1]
            _translate_floe!(FT, get_floe(floes, elem_idx), -Lx, zero(FT))
            _translate_floe!(FT, get_floe(floes, new_nfloes), Lx, zero(FT))
        end
    end
    return
end

#=
Find ghosts of given element and its known ghosts through a northern or
southern periodic boundary. If element's centroid isn't within the domain in the
north/south direction, swap it with its ghost since the ghost's centroid must
then be within the domain. `elem_idx`is the floe of floe of interest's index within
the floe list.
=#
function find_ghosts!(
    floes,
    elem_idx,
    nbound::PeriodicBoundary{North, FT},
    sbound::PeriodicBoundary{South, FT},
) where {FT <: AbstractFloat}
    Ly =  nbound.val - sbound.val
    nfloes = length(floes)
        # passing through southern boundary
    if (floes.centroid[elem_idx][2] - floes.rmax[elem_idx] < sbound.val)
        ghosts_on_bounds!(floes, elem_idx, sbound, [FT(0), Ly])
    # passing through northern boundary
    elseif (floes.centroid[elem_idx][2] + floes.rmax[elem_idx] > nbound.val) 
        ghosts_on_bounds!(floes, elem_idx, nbound, [FT(0), -Ly])
    end
    # if element centroid isn't in domain's north/south direction, swap with ghost
    new_nfloes = length(floes)
    if new_nfloes > nfloes
        if floes.centroid[elem_idx][2] < sbound.val
            _translate_floe!(FT, get_floe(floes, elem_idx), zero(FT), Ly)
            _translate_floe!(FT, get_floe(floes, new_nfloes), zero(FT), -Ly)
        elseif nbound.val < floes.centroid[elem_idx][2]
            _translate_floe!(FT, get_floe(floes, elem_idx), zero(FT), -Ly)
            _translate_floe!(FT, get_floe(floes, new_nfloes), zero(FT), Ly)
        end
    end
    return
end

#=
Add ghosts of all of the given floes passing through the two given boundaries to
the list of floes. The `max_boundary` is the northern or eastern boundary of domain.
The `min_boundary` is the southern or western boundary of domain.
=#
function add_floe_ghosts!(floes, max_boundary, min_boundary)
    nfloes = length(floes)
    # uses initial length of floes so we can append to list
    for i in eachindex(floes)
        # the floe is active in the simulation and a parent floe
        if floes.status[i].tag == active && floes.ghost_id[i] == 0
            find_ghosts!(
                floes,
                i,
                max_boundary,
                min_boundary,
            )
            new_nfloes = length(floes)
            if new_nfloes > nfloes
                nghosts = new_nfloes - nfloes
                # set ghost floe's ghost id
                floes.ghost_id[(nfloes + 1):new_nfloes] .= range(1, nghosts) .+ length(floes.ghosts[i])
                # remove ghost floes ghosts as these were added to parent
                empty!.(floes.ghosts[nfloes + 1:new_nfloes])
                # index of ghosts floes saved in parent
                append!(floes.ghosts[i], (nfloes + 1):(nfloes + nghosts))
                nfloes += nghosts
            end
        end
    end
    return
end

# When there are no periodic boundaries, no ghosts should be added.
function add_ghosts!(
    elems,
    ::Domain{
        FT,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
    },
) where {FT<:AbstractFloat}
    return
end

# Add ghosts for elements that pass through the northern or southern periodic boundaries. 
function add_ghosts!(
    elems,
    domain::Domain{
        FT,
        <:PeriodicBoundary,
        <:PeriodicBoundary,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
    },
) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.north, domain.south)
    return
end


# Add ghosts for elements that pass through the eastern or western periodic boundaries. 
function add_ghosts!(
    elems,
    domain::Domain{
        FT,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
        <:PeriodicBoundary,
        <:PeriodicBoundary,
    },
) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.east, domain.west)
    return
end

"""
    add_ghosts!(elems, domain)

Add ghosts for elements (floes of topography elements) that pass through any of the domain boundaries.
Add ghost elements to list in-place. Note that ghosts will only be added when the elements pass through
periodic boundaries.

## _Positional arguments_
- `elems::Union{StructArray{Floe}, StructArray{TopographyElement}`: list of elements (floes or topigraphy) to
    add "ghost" elements to.
- $DOMAIN_DEF
"""
function add_ghosts!(
    elems,
    domain::Domain{
        FT,
        <:PeriodicBoundary,
        <:PeriodicBoundary,
        <:PeriodicBoundary,
        <:PeriodicBoundary
    }
) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.east, domain.west)
    add_floe_ghosts!(elems, domain.north, domain.south)
    return
end

#=
Adding ghosts to floes on a device (`FixedWidthFloes`).

The ghosts are added after the floes in use (i.e. currently active floes), in the order
of their parents: first the ghosts through the east/west boundaries, then those through
the north/south boundaries, which also makes ghosts of the first ones.

Each pass:
1. `_count_ghosts!`: per floe, the number of ghosts it needs.
2. A scan of these counts, which gives each floe the rows of its new ghosts.
3. `_plan_ghosts!`: per floe, the source row, ghost id, and shift of each of its ghosts.
4. `_ghost_values_kernel!`: one thread per value of each ghost, which copies (and shifts)
   that value from the source row. One thread per floe would copy a few hundred values one
   after another, which is slow on a GPU, where there are often only a few hundred ghosts.
5. `_link_ghosts!`: per floe, add its ghosts to its list of ghosts, and swap places with its
   last ghost if its centroid is outside of the domain.
=#

#=
Scratch memory for adding ghosts to floes with the capacity of `floes`, on the same device:
the number of new ghosts per floe, its inclusive scan, AcceleratedKernels' scratch memory
for the scan, a host buffer for the total, and per new ghost its source row, ghost id, and
shift. Make it once and pass it to `add_ghosts!` so that it doesn't allocate device memory.
=#
struct GhostWorkspace{IT1, AT1, W}
    counts::IT1
    offsets::IT1
    scan::W
    total::Vector{Int32}
    sources::IT1
    ghost_ids::IT1
    shifts::AT1
end

function GhostWorkspace(floes::FixedWidthFloes)
    counts = similar(floes.num_ghosts)
    offsets = similar(counts)
    scan = AK.workspace(AK.accumulate!, +, offsets, counts)
    return GhostWorkspace(counts, offsets, scan, zeros(Int32, 1), similar(counts),
        similar(counts), similar(floes.rmax))
end

"""
    add_ghosts!(floes::FixedWidthFloes{FT}, domain::Domain{FT}, backend; workspace = nothing)

Add ghosts of the floes that pass through periodic boundaries of `domain` on `backend`.
The ghosts are added after the floes in use, and their indices are stored in the `ghosts`
of their parents. `domain` must have the float type of the floes, as the kernels don't
convert its values.

Returns the floes with the ghosts in use (see `set_n_floes`) and the workspace that fits
them (`nothing` if none was given and none was needed). `workspace` is a `GhostWorkspace`
for the capacity of `floes`, made if not given. This waits for the kernels to finish,
because the number of ghosts is needed on the host.

The returned floes share the arrays of `floes` unless they had to be made larger (see
`with_capacity`). Either way, replace `floes` with the returned floes, and pass the returned
workspace to the next call, so that it isn't made again. If the arrays are made larger in
the north/south pass, the arrays of `floes` are left half-updated, with only the east/west
ghosts.

The floes must not have ghosts yet (call `remove_ghosts!` first): `ghosts` has room for
`MAX_GHOSTS` ghosts per floe, which is the most a floe without ghosts can get, and the kernels
don't check it.
"""
function add_ghosts!(floes::FixedWidthFloes{FT}, domain::Domain{FT}, backend;
    workspace = nothing) where FT
    floes, workspace = add_floe_ghosts!(floes, domain.east, domain.west, backend, workspace)
    return add_floe_ghosts!(floes, domain.north, domain.south, backend, workspace)
end

# No ghosts through boundaries that aren't periodic
add_floe_ghosts!(floes::FixedWidthFloes{FT}, ::AbstractBoundary, ::AbstractBoundary, backend,
    workspace) where FT = floes, workspace

#=
Add the ghosts through `max_boundary` (north or east) and `min_boundary` (south or west).
Writes into `floes`, but returns the floes with the ghosts in use, which may have new,
larger arrays (see `with_capacity`), and the workspace that fits them.
=#
function add_floe_ghosts!(
    floes::FixedWidthFloes{FT},
    max_boundary::PeriodicBoundary,
    min_boundary::PeriodicBoundary,
    backend,
    workspace,
) where FT
    n = floes.n_floes
    n == 0 && return floes, workspace
    axis = _ghost_axis(max_boundary)
    min_val, max_val = min_boundary.val, max_boundary.val
    min_ext = GI.extent(min_boundary.poly)
    max_ext = GI.extent(max_boundary.poly)
    workspace = _ghost_workspace(workspace, floes)
    total = _count_ghosts(floes, axis, min_val, max_val, min_ext, max_ext, backend, workspace)
    total == 0 && return floes, workspace
    if n + total > capacity(floes)
        # rare: make room for the ghosts and 25% more, and count again with a new workspace
        floes = with_capacity(floes, n + total + cld(n + total, 4))
        workspace = GhostWorkspace(floes)
        _count_ghosts(floes, axis, min_val, max_val, min_ext, max_ext, backend, workspace)
    end
    L = max_val - min_val
    ws = workspace
    launch_per_floe!(_plan_ghosts!, backend, floes, ws.counts, ws.offsets, ws.sources,
        ws.ghost_ids, ws.shifts, axis, L, min_val)
    _ghost_values_kernel!(backend, (256, 1))(floes, n, ws.sources, ws.ghost_ids, ws.shifts,
        axis; ndrange = (total, _values_per_floe(floes)))
    launch_per_floe!(_link_ghosts!, backend, floes, ws.counts, ws.offsets, axis, L, min_val,
        max_val)
    KernelAbstractions.synchronize(backend)
    return set_n_floes(floes, n + total), workspace
end

#=
The coordinate that crosses the boundary: 1 (x) for east/west, 2 (y) for north/south.
It is passed to the kernels as a number.
=#
_ghost_axis(::PeriodicBoundary{<:Union{East, West}}) = 1
_ghost_axis(::PeriodicBoundary{<:Union{North, South}}) = 2

# The given workspace if it fits the capacity of floes, else a new one
function _ghost_workspace(workspace, floes)
    if isnothing(workspace) || length(workspace.counts) != capacity(floes)
        return GhostWorkspace(floes)
    end
    return workspace
end

#=
Count the ghosts each floe needs into `workspace.counts` and scan them into
`workspace.offsets`. Returns the total number of ghosts. The counts of the spare rows are
zero, so the last element of the scan is the total.
=#
function _count_ghosts(floes, axis, min_val, max_val, min_ext, max_ext, backend, workspace)
    fill!(workspace.counts, 0)
    launch_per_floe!(_count_ghosts!, backend, floes, workspace.counts, axis, min_val, max_val,
        min_ext, max_ext)
    # the scan may run on the host's threads, so the kernel must have finished
    KernelAbstractions.synchronize(backend)
    AK.accumulate!(+, workspace.offsets, workspace.counts; workspace = workspace.scan)
    last = length(workspace.offsets)
    copyto!(workspace.total, 1, workspace.offsets, last, 1)
    return Int(workspace.total[1])
end

#=
The number of ghosts that floe `i` needs: if an active parent floe overlaps the boundary
that it may cross, it needs a ghost of itself and of each of its ghosts.
Every vertex is within `rmax` of the centroid, so a floe can only overlap a boundary if its
centroid is within `rmax` of it. Therefore a quick test is performed to exclude some
polygons from the (more expensive) intersection test.

Note that a floe that is wider than the domain (i.e. crosses both boundaries) would not
be handled properly (although it would not be realistic, in any case).
=#
function _count_ghosts!(floes::FixedWidthFloes, i, counts, axis, min_val, max_val, min_ext, max_ext)
    (floes.status[i] == Int32(active) && floes.ghost_id[i] == 0) || return
    c = floes.centroid[i, axis]
    r = floes.rmax[i]
    if c - r < min_val
        GOG.interiors_intersect(_floe_poly(floes, i), min_ext) || return
    elseif c + r > max_val
        GOG.interiors_intersect(_floe_poly(floes, i), max_ext) || return
    else
        return
    end
    counts[i] = floes.num_ghosts[i] + Int32(1)
    return
end

#=
Plan the ghosts of floe `i` in the rows that the scan gave it: copies
of its ghosts (in order), then of itself, translated by the domain size `L` to the other
side of the domain, with ghost ids that follow those of floe i's ghosts. Ghost j (row
`n_floes + j`) is a copy of row `sources[j]`, shifted by `shifts[j]` along `axis`.
=#
function _plan_ghosts!(floes::FixedWidthFloes, i, counts, offsets, sources, ghost_ids, shifts, axis, L, min_val)
    nnew = counts[i]
    nnew == 0 && return
    before = offsets[i] - nnew
    # through the min boundary, the ghost is on the max side of the domain, and vice versa
    shift = floes.centroid[i, axis] - floes.rmax[i] < min_val ? L : -L
    nghosts = floes.num_ghosts[i]
    # Int32, as the index arrays: the stores don't convert
    for k in Int32(1):nnew
        j = before + k
        sources[j] = k <= nghosts ? floes.ghosts[i, k] : Int32(i)
        ghost_ids[j] = nghosts + k
        shifts[j] = shift
    end
    return
end

#=
Write value v of ghost j: one thread per value, so that a ghost's values are written at the
same time, not one after another.
=#
@kernel function _ghost_values_kernel!(floes, n, sources, ghost_ids, shifts, axis)
    j, v = @index(Global, NTuple)
    @inline _ghost_value!(floes, n + j, sources[j], v, axis, shifts[j], ghost_ids[j])
end

#=
Fields that a ghost does not copy from its parent: a new ghost has no interactions,
no damage (?) and no ghosts of its own, and flags are for the events of the floe itself.

TODO: Check if it is correct for a new ghost to have no damage - this seems to be how
the original CPU version of Subzero did it.
=#
const _GHOST_ZEROED_FIELDS = (:num_inters, :interactions, :damage, :num_ghosts, :ghosts, :flags)

#=
Value v of row dst, the ghost of row src: the values of a row are numbered field by field,
in the order of the fields of `FixedWidthFloes`, and within a field in memory order. The
value is copied from src, and shifted by `shift` along `axis` for the centroid and the
vertices of the polygon (`src + shift`, as `_translate_floe!`), except for the fields in
`_GHOST_ZEROED_FIELDS`, which are zero, and the ghost id. Generated, so that the field and
what to do with it are known at compile time.
=#
@generated function _ghost_value!(floes::FixedWidthFloes, dst, src, v, axis, shift, ghost_id)
    body = Expr(:block, Expr(:meta, :inline), :(first = 0))
    for (k, name) in enumerate(fieldnames(floes))
        fieldtype(floes, k) <: AbstractArray || continue
        value = if name in _GHOST_ZEROED_FIELDS
            :(zero(eltype(a)))
        elseif name == :ghost_id
            :(ghost_id)
        elseif name == :centroid
            :(_element(a, src, e) + (e == axis ? shift : zero(shift)))
        elseif name == :poly
            # e is vertex p, coordinate d; the padding after the last vertex stays zero
            quote
                p, d = mod1(e, size(a, 2)), cld(e, size(a, 2))
                p <= floes.n_points[src] ?
                    _element(a, src, e) + (d == axis ? shift : zero(shift)) :
                    _element(a, src, e)
            end
        else
            :(_element(a, src, e))
        end
        push!(body.args, quote
            a = getfield(floes, $k)
            len = _row_length(a)
            if first < v <= first + len
                e = v - first
                _set_element!(a, dst, e, $value)
                return
            end
            first += len
        end)
    end
    push!(body.args, :(return))
    return body
end

# Number of values in a row of an array of a FixedWidthFloes, and element e of row i, with
# the elements of a row in memory order. Inlined, so that the arrays aren't passed by
# reference, which made CUDA copy all of the floes into each thread's local memory.
_row_length(a::AbstractArray{<:Any, 1}) = 1
_row_length(a::AbstractArray{<:Any, 2}) = size(a, 2)
_row_length(a::AbstractArray{<:Any, 3}) = size(a, 2) * size(a, 3)
@inline _element(a::AbstractArray{<:Any, 1}, i, e) = a[i]
@inline _element(a::AbstractArray{<:Any, 2}, i, e) = a[i, e]
@inline _element(a::AbstractArray{<:Any, 3}, i, e) =
    a[i, mod1(e, size(a, 2)), cld(e, size(a, 2))]
@inline _set_element!(a::AbstractArray{<:Any, 1}, i, e, x) = (a[i] = x; nothing)
@inline _set_element!(a::AbstractArray{<:Any, 2}, i, e, x) = (a[i, e] = x; nothing)
@inline _set_element!(a::AbstractArray{<:Any, 3}, i, e, x) =
    (a[i, mod1(e, size(a, 2)), cld(e, size(a, 2))] = x; nothing)

# Number of values in a row of floes, over all fields
_values_per_floe(floes::FixedWidthFloes) = sum(fieldnames(FixedWidthFloes)) do name
    field = getfield(floes, name)
    field isa AbstractArray ? _row_length(field) : 0
end

#=
Add the new ghosts of floe i to its list of ghosts. If floe i's centroid is outside of the
domain, swap places with its last new ghost, the copy of floe i, so that the floe in use is
the one inside the domain. Note that this "swap" is a swap of coordinates, not an exchange
of rows in the Floes array.
=#
function _link_ghosts!(floes::FixedWidthFloes, i, counts, offsets, axis, L, min_val, max_val)
    nnew = counts[i]
    nnew == 0 && return
    # rows of the new ghosts: before + 1, ..., before + nnew, as Int32 like the index arrays
    before = Int32(floes.n_floes) + offsets[i] - nnew
    nghosts = floes.num_ghosts[i]
    for k in Int32(1):nnew
        floes.ghosts[i, nghosts + k] = before + k
    end
    floes.num_ghosts[i] = nghosts + nnew
    c = floes.centroid[i, axis]
    if c < min_val
        _translate_floe!(floes, i, _along_axis(axis, L)...)
        _translate_floe!(floes, before + nnew, _along_axis(axis, -L)...)
    elseif max_val < c
        _translate_floe!(floes, i, _along_axis(axis, -L)...)
        _translate_floe!(floes, before + nnew, _along_axis(axis, L)...)
    end
    return
end

# A translation by Δ along an axis (1 = x, 2 = y), as (Δx, Δy)
_along_axis(axis, Δ) = axis == 1 ? (Δ, zero(Δ)) : (zero(Δ), Δ)

#=
Remove the ghosts that `add_ghosts!` added: keep the first `n_init_floes` floes, which
are the floes without ghosts, and empty their lists of ghosts. Floes that
swapped places with a ghost keep their new place.

Returns the floes without ghosts, which share the arrays of `floes`: replace `floes` with
them. The kernel is launched asynchronously. The rows of the ghosts keep their values (they
aren't zero and their status isn't 0), so fill them before using them.
=#
function remove_ghosts!(floes::FixedWidthFloes, n_init_floes, backend)
    floes = set_n_floes(floes, n_init_floes)
    launch_per_floe!(_clear_ghosts!, backend, floes)
    return floes
end

function _clear_ghosts!(floes::FixedWidthFloes, i)
    floes.num_ghosts[i] = Int32(0)
    _zero_row!(floes.ghosts, i)
    return
end
