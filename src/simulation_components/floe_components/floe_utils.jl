const FLOES_DEF = "`floes::StructArray{Floe}`: simulation floes"

# Syntactic sugar for use in code
const FloeType{FT} = Union{LazyRow{Floe{FT}}, Floe{FT}} where FT

# lazily access floe within struct array floe field
get_floe(floes::StructArray, i::Int) = LazyRow(floes, i)

#=
Split a given polygon around any holes before turning each region with an area greater than
the minimum floe area into a floe. New floes are added to existing loe field.
Any Floe constructed keyword arguments can be passed after required arguments.
Returns the total number of floes created.
=#
function _poly_to_floes!(::Type{FT}, floes, poly, hmean, Δh, rmax;
    floe_settings, rng = Xoshiro(), kwargs...
) where {FT <: AbstractFloat}
    a = area_poly(poly, FT)
    # only make polygon into floe if it is big enough
    if a >= floe_settings.min_floe_area && a > 0
        # if it doesn't have a hole, add right to list
        if !hashole(poly)
             height = clamp( # set floe height
                hmean + (-1)^rand(rng, 0:1) * rand(rng, FT) * Δh,
                floe_settings.min_floe_height,
                floe_settings.max_floe_height,
            )
            floe = Floe{FT}( # create floe
                poly::Polys,
                height;
                floe_settings,
                rng,
                kwargs...
            )
            push!(floes, floe) # add to floe field
            return 1 # one floe added
        else
            # split floe around first hole
            cx, cy = centroid_poly(GI.gethole(poly, 1), FT)
            new_regions = cut_poly_by_line(poly, GI.Line([(cx - rmax, cy), (cx + rmax, cy)]), FT)
            n = 0
            for r in new_regions # recrusively call function around two new pieces
                n += _poly_to_floes!(FT, floes, r, hmean, Δh, rmax;
                    floe_settings, rng, kwargs...)
            end
            return n
        end
    end
    return 0
end

# translate floe by Δx and Δy
function _translate_floe!(::Type{FT}, floe, Δx, Δy) where FT
    floe.centroid[1] += Δx
    floe.centroid[2] += Δy
    floe.poly = _translate_poly(FT, floe.poly, Δx, Δy)
    floe.coords = find_poly_coords(floe.poly, FT)
    return
end

# translate floe by floe by Δx and Δy and rotate flow by Δα
function _move_floe!(::Type{FT}, floe, Δx, Δy, Δα) where FT
    cx, cy = floe.centroid
    # move centroid
    floe.centroid[1] += Δx
    floe.centroid[2] += Δy
    # move polygon
    floe.poly = _move_poly(FT, floe.poly, Δx, Δy, Δα, cx, cy)::Polys{FT}
    floe.coords = find_poly_coords(floe.poly, FT)
    return 
end


#= 
Deepcopy of a floe by creating a new floe and deep copying all fields.
New floe with fields that are equal in value. Any vector fields are copies so
they share values, but not memory location.
=#
function deepcopy_floe(floe::LazyRow{Floe{FT}}) where {FT}
    poly = get_tuple_poly(floe.poly, FT)
    f = Floe{FT}(
        poly = poly,
        centroid = copy(floe.centroid),
        coords = find_poly_coords(poly, FT),
        height = floe.height,
        area = floe.area,
        mass = floe.mass,
        rmax = floe.rmax,
        moment = floe.moment,
        angles = copy(floe.angles),
        bonds = copy(floe.bonds),
        bond_break_mode = copy(floe.bond_break_mode),
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
        fxA = floe.fxA,
        fyA = floe.fyA,
        fxO = floe.fxO,
        fyO = floe.fyO,
        fxP = floe.fxP,
        fyP = floe.fyP,
        fxC = floe.fxC,
        fyC = floe.fyC,
        trqOA = floe.trqOA,
        hflx_factor = floe.hflx_factor,
        overarea = floe.overarea,
        collision_force = copy(floe.collision_force),
        collision_trq = floe.collision_trq,
        bond_force = copy(floe.bond_force),
        bond_trq = floe.bond_trq,
        stress_bond = copy(floe.stress_bond),
        ind_bond_force_nt = deepcopy(floe.ind_bond_force_nt),
        ind_bond_force_xy = deepcopy(floe.ind_bond_force_xy),
        ind_bond_stiff_nt = deepcopy(floe.ind_bond_stiff_nt),
        ind_bond_normal_old = deepcopy(floe.ind_bond_normal_old),
        ind_bond_bm = copy(floe.ind_bond_bm),
        ind_bond_prop = deepcopy(floe.ind_bond_prop),
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

#=
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
    xc, yc = get_tuple_point(cent, T)
    Ixx, Iyy = zero(T), zero(T)
    x1, y1 = zero(T), zero(T)
    for (i, p2) in enumerate(GI.getpoint(GI.getexterior(poly)))
        (x2, y2) = get_tuple_point(p2, T)
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

# ensure provided floe coordinates (of type PolyVec) are valid and without any holes
function _correct_floe_poly(::Type{FT}, coords::PolyVec) where FT
    valid_polyvec!(coords)
    poly = make_polygon(coords, FT)
    return _correct_floe_poly(FT, poly)
end
function signed_area_ring(ring)   # ring: vector of (x, y) or [x, y]
    s = 0.0
    for i in 1:length(ring)-1     # closing point repeats the first, so stop one early
        s += ring[i][1] * ring[i+1][2] - ring[i+1][1] * ring[i][2]
    end
    return s / 2
end
# ensure provided polygon points are of the right type and that the polygon has no holes
function _correct_floe_poly(::Type{FT}, poly::Polys) where FT
    poly = get_tuple_poly(poly, FT)
    rmholes!(poly)
    # enforce CCW on the exterior ring
    coords = find_poly_coords(poly, FT)
    if signed_area_ring(coords[1]) < 0
        reverse!(coords[1])       # closed ring: reversing the whole thing keeps first == last
        poly = make_polygon(coords, FT)
    end
    return poly
end

struct BoundingBox
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end

function compute_bbox(poly_coords::Vector{Vector{Float64}})
    xs = [c[1] for c in poly_coords]
    ys = [c[2] for c in poly_coords]
    BoundingBox(minimum(xs), maximum(xs), minimum(ys), maximum(ys))
end

function bbox_intersect(b1::BoundingBox, b2::BoundingBox; atol=1e-8)
    return !(
        (b1.xmax <= b2.xmin - atol) ||  (b1.xmin >= b2.xmax + atol) ||
        (b1.ymax <= b2.ymin - atol) ||  (b1.ymin >= b2.ymax + atol))
end

# Do two edges (a→b on floe i, c→d on floe j) coincide within tolerance?
# With both floes CCW, a shared edge is traversed in opposite directions: a≈d and b≈c.
@inline function _edges_match(a, b, c, d, rtol)
    tol = rtol * 0.5 * (hypot(b[1] - a[1], b[2] - a[2]) + hypot(d[1] - c[1], d[2] - c[2]))
    return abs(a[1] - d[1]) ≤ tol && abs(a[2] - d[2]) ≤ tol &&
           abs(b[1] - c[1]) ≤ tol && abs(b[2] - c[2]) ≤ tol
end

# Bond every pair of matching edges between floes i and j
function _bond_shared_edges!(floe_arr, i, j, ring_i, ring_j, rtol)
    bonds_i, bonds_j = floe_arr.bonds[i], floe_arr.bonds[j]
    for ei in 1:length(ring_i) - 1
        bonds_i[ei] != 0 && continue
        a, b = ring_i[ei], ring_i[ei + 1]
        for ej in 1:length(ring_j) - 1
            bonds_j[ej] != 0 && continue
            if _edges_match(a, b, ring_j[ej], ring_j[ej + 1], rtol)
                bonds_i[ei] = floe_arr.id[j]
                bonds_j[ej] = floe_arr.id[i]
                break
            end
        end
    end
    return
end

function calculate_bonds!(floe_arr; rtol = 1e-3, bbox_atol = 1e-8)
    n = length(floe_arr)
    rings  = [floe_arr.coords[i][1] for i in 1:n]   # exterior rings, built once
    bboxes = [compute_bbox(r) for r in rings]       # built once
    foreach(b -> fill!(b, 0), floe_arr.bonds)       # start clean
    # Sweep along x: only compare floes whose x-ranges overlap
    order = sortperm(bboxes; by = b -> b.xmin)
    for a in 1:n
        i = order[a]
        for b in a + 1:n
            j = order[b]
            bboxes[j].xmin > bboxes[i].xmax + bbox_atol && break   # all later floes are further right
            bbox_intersect(bboxes[i], bboxes[j]; atol = bbox_atol) || continue
            _bond_shared_edges!(floe_arr, i, j, rings[i], rings[j], rtol)
        end
    end
    return floe_arr
end
