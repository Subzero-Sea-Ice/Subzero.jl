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
    return
end

# translate floe by floe by Δx and Δy and rotate flow by Δα
function _move_floe!(floes, i, Δx, Δy, Δα) 
    cx, cy = floes.centroid[i, :]
    # move centroid
    floes.centroid[i, 1] += Δx
    floes.centroid[i, 2] += Δy
    # move polygon
    n_points = 0
    for x in floes.poly[i, :, 1]
        if x == FILL_VALUE
            break
        end
        n_points += 1
    end
    points = [Tuple(floes.poly[i, j, :]) for j in 1:n_points]
    polygon = GI.Polygon([GI.LinearRing(points)])
    moved_polygon = _move_poly(polygon, Δx, Δy, Δα, cx, cy)
    moved_points = GI.getpoint(moved_polygon)
    for (j, point) in enumerate(moved_points)
        floes.poly[i, j, :] .= point
    end
    points_dim = 2
    for j in n_points+1:size(floes.poly, points_dim)
        floes.poly[i, j, :] .= FILL_VALUE
    end
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

# ensure provided polygon points are of the right type and that the polygon has no holes
function _correct_floe_poly(::Type{FT}, poly::Polys) where FT
    poly = get_tuple_poly(poly, FT)
    rmholes!(poly)
    return poly
end