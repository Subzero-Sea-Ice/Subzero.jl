# Functions needed for fracturing floes

"""
    determine_fractures(...)

Determines which floes will fracture depending on the principal stress criteria.
Runs [`determine_fractures!`](@ref) on the CPU.

## _Positional arguments_
- $FLOES_DEF
- `criteria::AbstractFractureCriteria`: floe fracture criteria
- $FLOE_SETTINGS_DEF

##  _Returns_
- `Vector{Int}`: list of indices of floes to fracture
"""
function determine_fractures(
    floes::StructArray{<:Floe},
    criteria::AbstractFractureCriteria,
    floe_settings,
)
    fwf = FixedWidthFloes(floes)
    determine_fractures!(fwf, criteria, floe_settings; backend = CPU())
    return _flagged(fwf.flags, FLAG_FRACTURE)
end

"""
    determine_fractures!(floes, criteria, floe_settings; backend = CPU())

Determines which floes will fracture depending on the principal stress criteria, and sets
`FLAG_FRACTURE` in `floes.flags` for those floes. A floe fractures if its accumulated stress
in principal stress space (see [`find_σpoint`](@ref)) is outside of the criteria polygon and
its area is at least `floe_settings.min_floe_area`.

The floes must already be on `backend`. The kernel is launched asynchronously, so call
`KernelAbstractions.synchronize(backend)` before reading the flags.

## _Positional arguments_
- `floes::FixedWidthFloes`: floes on `backend`
- `criteria::AbstractFractureCriteria`: floe fracture criteria, updated in-place with
    `_update_criteria!(criteria, floes)`
- $FLOE_SETTINGS_DEF

## _Keyword arguments_
- `backend::KernelAbstractions.Backend`: backend to run the calculations on (Default = `CPU()`)
"""
function determine_fractures!(
    floes::FixedWidthFloes{FT},
    criteria::AbstractFractureCriteria,
    floe_settings;
    backend = CPU(),
) where FT
    #= Update the criteria with the floes' current values (e.g. their mean height). Reductions
    such as `mean` run after the kernels already launched on the backend, so no synchronize
    is needed. =#
    _update_criteria!(criteria, floes)
    criteria_coords = adapt(backend, _criteria_coords(criteria.poly, FT))
    launch_per_floe!(determine_fracture!, backend, floes, criteria_coords,
        FT(floe_settings.min_floe_area), floe_settings.stress_calculator)
    return
end

# Coordinates of the exterior of the criteria polygon in principal stress space, as a
# (number of points, 2) matrix that can be moved to a device.
function _criteria_coords(poly, ::Type{FT}) where FT
    ring = GI.getexterior(poly)
    coords = Matrix{FT}(undef, GI.npoint(ring), 2)
    for (j, point) in enumerate(GI.getpoint(ring))
        coords[j, 1], coords[j, 2] = get_tuple_point(point, FT)
    end
    return coords
end

#=
Sets `FLAG_FRACTURE` for floe i of a `FixedWidthFloes` if it should fracture: its stress
point is outside of the criteria polygon `criteria_coords` and its area is at least
`min_floe_area`. Clears the flag otherwise.
=#
function determine_fracture!(floes::FixedWidthFloes, i, criteria_coords, min_floe_area, stress_calculator)
    σ1, σ2 = find_σpoint(floes, i, stress_calculator, min_floe_area)
    fractures = floes.area[i] >= min_floe_area &&
        !_point_covered_by(σ1, σ2, criteria_coords, size(criteria_coords, 1))
    floes.flags[i] = fractures ? floes.flags[i] | FLAG_FRACTURE : floes.flags[i] & ~FLAG_FRACTURE
    return
end

"""
    find_σpoint(floes::FixedWidthFloes, i, stress_calculator, min_floe_area)

Find floe `i`'s accumulated stress in principal stress space, so that it can be compared to
the fracture criteria: the eigenvalues of `stress_accum` in ascending order, scaled with
`_scale_principal_stress`. Written as scalar operations, so it can be called from a GPU
kernel.

##  _Returns_
- `Tuple{FT, FT}`: the two scaled principal stresses
"""
function find_σpoint(floes::FixedWidthFloes{FT}, i, stress_calculator, min_floe_area) where FT
    a = floes.stress_accum[i, 1, 1]
    b = floes.stress_accum[i, 1, 2]
    c = floes.stress_accum[i, 2, 1]
    d = floes.stress_accum[i, 2, 2]
    #= Eigenvalues of [a b; c d], which calc_stress! keeps symmetric (b == c). hypot doesn't
    underflow for tiny stresses, unlike sqrt(((a - d) / 2)^2 + b^2), which would make both
    eigenvalues equal. =#
    mid = (a + d) / 2
    r = hypot((a - d) / 2, (b + c) / 2)
    return _scale_principal_stress(stress_calculator, mid - r, mid + r, floes, i, min_floe_area)
end

"""
    deform_floe!(...)

Deform a floe around the area of its collision with largest area overlap within
the last timestep.

The input floe's centroid, coordinates, and area are updated to reflect a deformation
due to the collision with the deforming floe.

## _Positional arguments_
- $FLOE_DEF
- `deformer_poly::Polys`: polygon of floe that is deforming floe argument
- `deforming_forces:: Vector{Float}`: 1x2 matrix of forces between floe and
    the deforming floe from floe's interactions - of the form: [xforce yforce]
- $FLOE_SETTINGS_DEF
- $ΔT_DEF
- `rng::RandomNumberGenerator`:: random number generator
"""
function deform_floe!(
    floe,
    deformer_poly::Polys{FT},
    deforming_forces,
    floe_settings,
    Δt,
    rng,
) where FT
    poly = floe.poly
    overlap_regions = intersect_polys(poly, deformer_poly, FT)
    FT_area_poly(p) = area_poly(p, FT)
    max_overlap_area, max_overlap_idx = findmax(FT_area_poly, overlap_regions)
    overlap_region = overlap_regions[max_overlap_idx]
    # If floe and the deformer floe have an overlap area
    if max_overlap_area > 0
        # Determine displacement of deformer floe
        region_cent = centroid_poly(overlap_region, FT)
        dist = dist_to_poly(region_cent,overlap_region, FT)
        force_fracs = deforming_forces ./ 2norm(deforming_forces)
        Δx, Δy = abs.(dist)[1] .* force_fracs
        # Temporarily move deformer floe to find new shape of floe
        deformer_poly = _translate_poly(FT, deformer_poly, Δx, Δy)
        new_floes = diff_polys(poly, deformer_poly, FT)
        new_floe_area, new_floe_idx = findmax(FT_area_poly, new_floes)
        new_floe_poly = new_floes[new_floe_idx]
        # If didn't change floe area by more than 90%
        if new_floe_area > 0 && new_floe_area/floe.area > 0.9
            # Update floe shape and conserve mass
            moment_tmp = floe.moment
            x_tmp, y_tmp = floe.centroid
            replace_floe!(
                floe,
                new_floe_poly,
                floe.mass,
                floe_settings,
                rng,
            )
            conserve_momentum_change_floe_shape!(
                floe.mass,
                moment_tmp,
                x_tmp,
                y_tmp,
                Δt,
                floe,
            )
        end
    end
    return
end

"""
    split_floe(...)

Splits a given floe into pieces using voronoi tesselation.
User will recieve a warning if floe isn't split.

## _Positional arguments_
- $FLOE_DEF
- `rng::RandomNumberGenerator`:: random number generator
- `fracture_settings::FractureSettings`: simulation's fracture settings
- $FLOE_SETTINGS_DEF
- $ΔT_DEF
 
##  _Returns_
- `new_floes::StructArray{Floes}`: list of pieces floe is split into, each of which is a new floe
"""
function split_floe(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    rng,
    fracture_settings,
    floe_settings,
    Δt,
) where {FT}
    new_floes = StructArray{Floe{FT}}(undef, 0)
    # Generate voronoi tesselation in floe's bounding box
    xmin, ymin = floe.centroid[1] - floe.rmax, floe.centroid[2] - floe.rmax
    Δx, Δy = 2floe.rmax, 2floe.rmax
    pieces_polys = _generate_voronoi_coords(FT,
        fracture_settings.npieces,
        Δx, Δy,
        xmin, ymin,
        floe.poly,
        rng,
        1,  # Warn if only 1 point is identified as the floe won't be split
    )
    if !isempty(pieces_polys)
        # Intersect voronoi tesselation pieces with floe
        rmholes!(floe.poly)
        pieces_polys = mapreduce(p -> intersect_polys(p, floe.poly, FT), append!, pieces_polys; init = Vector{Polys{FT}}())
        # Conserve mass within pieces
        pieces_areas = [area_poly(p, FT) for p in pieces_polys]
        total_area = sum(pieces_areas)
        # Create floes out of each piece
        for i in eachindex(pieces_polys)
            if pieces_areas[i] > 0
                mass = floe.mass * (pieces_areas[i]/total_area)
                height = mass / (floe_settings.ρi * pieces_areas[i])
                _poly_to_floes!(
                    FT,
                    new_floes,
                    pieces_polys[i],
                    height,
                    0, # Δh - range of random height difference between floes
                    floe.rmax;
                    floe_settings = floe_settings,
                    rng = rng,
                    u = floe.u,
                    v = floe.v,
                    ξ = floe.ξ,
                )
            end
        end
    end

    # Conserve momentum and update strain
    conserve_momentum_fracture_floe!(
        floe,
        new_floes,
        Δt,
    )
    for s in new_floes.strain
        s .= floe.strain
    end

    return new_floes
end

"""
    fracture_floes!(...)

Fractures floes that meet the criteria defined in the fracture settings.
Floe pieces added to floe array and original fractured floes are removed.
The new highest floe ID after adding the new floes to the array is returned.

## _Positional arguments_
- $FLOES_DEF
- `max_floe_id::Int`: maximum ID of any floe created so far in simulation
- `rng::RandomNumberGenerator`:: random number generator
- `fracture_settings::FractureSettings`: simulation's fracture settings
- $FLOE_SETTINGS_DEF
- $ΔT_DEF

## _Keyword arguments_
- `frac_idx::Vector{Int}`: indices of the floes to fracture, e.g. from the flags set by
    [`determine_fractures!`](@ref) (Default: calculated with [`determine_fractures`](@ref))

##  _Returns_
- `max_floe_id::Int`: new highest floe ID after adding new floes to floe array.
"""
function fracture_floes!(
    floes::StructArray{<:Floe{FT}},
    max_floe_id,
    rng,
    fracture_settings,
    floe_settings,
    Δt;
    frac_idx = determine_fractures(floes, fracture_settings.criteria, floe_settings),
) where {FT <: AbstractFloat}
    # Initialize list for new floes created from fracturing existing floes
    nfloes2frac = length(frac_idx)
    fracture_list = [StructArray{Floe{FT}}(undef, 0) for _ in 1:nfloes2frac]
    # Fracture floes that meet criteria 
    Threads.@threads for i in 1:nfloes2frac
        # Deform floe around largest impact site
        if fracture_settings.deform_on
            fidx = frac_idx[i]
            max_overlap = FT(0)
            max_overlap_idx = 0
            for i in 1:floes.num_inters[fidx]
                if (
                    floes.interactions[fidx][i, floeidx] > 0 &&
                    max_overlap < floes.interactions[fidx][i, overlap]
                )
                    max_overlap = floes.interactions[fidx][i, overlap]
                    max_overlap_idx = i
                end
            end
            if max_overlap > 0
                deforming_inter = floes.interactions[fidx][max_overlap_idx, :]
                deforming_floe_idx = Int(deforming_inter[floeidx])
                if deforming_floe_idx <= length(floes)
                    deform_floe!(
                        get_floe(floes, frac_idx[i]), 
                        floes.poly[deforming_floe_idx],
                        deforming_inter[xforce:yforce],
                        floe_settings,
                        Δt,
                        rng,
                    )
                end
            end
        end
        # Split flie into pieces
        new_floes = split_floe(
            get_floe(floes, frac_idx[i]),
            rng,
            fracture_settings,
            floe_settings,
            Δt,
        )
        append!(fracture_list[i], new_floes)
    end
    # Remove old (unfractured) floes and add fractured pieces
    for i in range(nfloes2frac, 1, step = -1)
        new_floes = fracture_list[i]
        if !isempty(new_floes)
            n_new_floes = length(new_floes)
            new_floes.id .= range(max_floe_id + 1, max_floe_id + n_new_floes)
            push!.(new_floes.parent_ids, floes.id[frac_idx[i]])
            append!(floes, new_floes)
            max_floe_id += n_new_floes
            StructArrays.foreachfield(f -> deleteat!(f, frac_idx[i]), floes)
        end
    end
    return max_floe_id
end
