"""
Functions to simplify and reduce individual floes and floe list. 
"""

"""
    dissolve_floe(floe, grid, dissolved)

Dissolve given floe into dissolved ocean matrix.
Inputs:
    floe        <Union{Floe, LazyRow{Floe}}> single floe
    grid        <RegRectilinearGrid> model's grid
    domain      <Domain> model's domain
    dissolved    <Matrix{AbstractFloat}> ocean's dissolved field
Outputs:
    None. Update dissolved matrix with given floe's mass and mark floe for
    removal.
"""
function dissolve_floe!(floe, grid::RegRectilinearGrid, domain, dissolved)
    xidx, yidx = find_grid_cell_index(
        floe.centroid[1],
        floe.centroid[2],
        grid,
    )
    # Wrap indexes around grid if bounds are periodic 
    xidx = shift_cell_idx(xidx, grid.Nx + 1, domain.east)
    yidx = shift_cell_idx(yidx, grid.Ny + 1, domain.north)
    # If centroid is within bounds after wrapping, add mass to dissolved
    if 0 < xidx <= grid.Nx && 0 < yidx <= grid.Ny
        dissolved[yidx, xidx] += floe.mass
    end
    return
end

"""
    smooth_floes!(
        floes,
        topography,
        simp_settings,
        collision_settings,
        Δt,
        rng,
    )
Smooths floe coordinates for floes with more vertices than the maximum
allowed number. Uses Ramer–Douglas–Peucker algorithm with a user-defined
tolerance. If new shape causes overlap greater with another floe greater than
the maximum percentage allowed, mark the two floes for fusion.
Inputs:
    floes               <StructArray{Floe}> model's floes
    topography          <StructArray{TopographyElement}> domain's topography
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    collision_settings  <CollisionSettings> simulation's collision settings
    Δt                  <Int> length of simulation timestep in seconds
    rng                 <RNG> random number generator for new monte carlo points
"""
function smooth_floes!(
    floes::StructArray{Floe{FT}},
    topography,
    simp_settings,
    collision_settings,
    floe_settings,
    Δt,
    rng,
) where {FT <: AbstractFloat}
    for i in eachindex(floes)
        if GI.npoint(GI.getexterior(floes.poly[i])) > simp_settings.max_vertices
            poly_list = [simplify_poly(floes.poly[i], simp_settings.tol)]
            if !isempty(topography)
                poly_list = diff_polys(make_multipolygon(poly_list), make_multipolygon(topography.poly), FT)
            end
            simp_poly =
                if length(poly_list) == 1
                    poly_list[1]
                else
                    areas = [GO.area(p) for p in poly_list]
                    _, max_idx = findmax(areas)
                    poly_list[max_idx]
                end
            x_tmp, y_tmp = floes.centroid[i]
            moment_tmp = floes.moment[i]
            replace_floe!(
                get_floe(floes, i),
                simp_poly,
                floes.mass[i],
                floe_settings,
                rng,
            )
            # conserve momentum
            conserve_momentum_change_floe_shape!(
                floes.mass[i],
                moment_tmp,
                x_tmp,
                y_tmp,
                Δt,
                get_floe(floes, i),
            )
            # Mark interactions for fusion
            for j in eachindex(floes)
                if i != j && floes.status[j].tag != remove && potential_interaction(
                    floes.centroid[i],
                    floes.centroid[j],
                    floes.rmax[i],
                    floes.rmax[j],
                )
                    if floes.status[j].tag == fuse && i in floes.status[j].fuse_idx
                        floes.status[i].tag = fuse
                        push!(floes.status[i].fuse_idx, j)
                    else
                        jpoly = floes.poly[j]
                        intersect_area = sum(GO.area, intersect_polys(simp_poly, jpoly, FT); init = 0.0)
                        if intersect_area/GO.area(jpoly) > collision_settings.floe_floe_max_overlap
                            floes.status[i].tag = fuse
                            push!(floes.status[i].fuse_idx, j)
                        end
                    end
                end
            end
        end
    end
    return
end

"""
    fuse_two_floes!(
        keep_floe,
        remove_floe,
        Δt,
        floe_settings,
        max_floe_id,
        rng,
    )
Fuses two floes together if they intersect and replaces the larger of the two
floes with their union. Mass and momentum are conserved.
Inputs:
    keep_floe           <Union{Floe, LazyRow{Floe}}> first floe
    remove_floe         <Union{Floe, LazyRow{Floe}}> second floe
    floe_settings       <FloeSettings> simulation's settings to make new floes
    prefuse_max_floe_id <Int> maximum floe ID used yet in simulation
    rng                 <RNG> random number generator
Outputs:
    If floes are not intersecting, no changes. If intersecing, the fused floe
    replaces the larger of the two floes and the smaller floe is marked for
    removal.
    Note that the smaller floe's ID is NOT updated!
"""
function fuse_two_floes!(
    keep_floe::FloeType{FT},
    remove_floe,
    Δt,
    floe_settings,
    prefuse_max_floe_id,
    rng,
) where FT
    # Create new polygon if they fuse
    rmholes!(keep_floe.poly)
    rmholes!(remove_floe.poly)
    poly1 = keep_floe.poly
    poly2 = remove_floe.poly
    new_poly_list = union_polys(poly1, poly2, FT)
    if length(new_poly_list) == 1  # if they fused, they will make one polygon
        new_poly = new_poly_list[1]
        rmholes!(new_poly)
        # mark smaller floe for removal
        remove_floe.status.tag = remove
        # record as value will change with replace
        mass_tmp = keep_floe.mass
        moment_tmp = keep_floe.moment
        x_tmp, y_tmp = keep_floe.centroid
        # give keep_floe new shape
        replace_floe!(
            keep_floe,
            new_poly,
            keep_floe.mass + remove_floe.mass,
            floe_settings,
            rng,
        )
        # conserve momentum
        conserve_momentum_change_floe_shape!(
            mass_tmp,
            moment_tmp,
            x_tmp,
            y_tmp,
            Δt,
            keep_floe,
            remove_floe,
        )
        # Update stress history
        keep_floe.stress_accum .= (1/keep_floe.mass) * (
            keep_floe.stress_accum * mass_tmp .+
            remove_floe.stress_accum * remove_floe.mass
        )
        keep_floe.stress_instant = (1/keep_floe.mass) * (
            keep_floe.stress_instant * mass_tmp + 
            remove_floe.stress_instant * remove_floe.mass)
        # Update IDs
        if 0 < keep_floe.id <= prefuse_max_floe_id
            push!(keep_floe.parent_ids, keep_floe.id)
        end
        if 0 < remove_floe.id <= prefuse_max_floe_id
            push!(keep_floe.parent_ids, remove_floe.id)
        end
    end
    return
end

"""
    fuse_floes!(
        floes,
        max_floe_id,
        floe_settings,
        Δt,
        rng,
    )

Fuse all floes marked for fusion.
Inputs:
    floes               <StructArray{Floe}> model's floes
    max_floe_id         <Int> maximum floe ID created yet
    floe_settings       <FloeSettings>  simulation's settings for making floes
    Δt                  <Int> simulation timestep in seconds
    rng                 <RNG> random number generator
Outputs:
    None. Fuses floes marked for fusion. Marks floes fused into another floe
    for removal. 
"""
function fuse_floes!(
    floes,
    max_floe_id,
    floe_settings,
    Δt,
    rng,
)
    prefuse_max_floe_id = max_floe_id
    for i in eachindex(floes)
        if floes.status[i].tag == fuse
            for j in floes.status[i].fuse_idx
                # not already fused with another floe or marked for removal
                if floes.status[j].tag != remove
                    keep_idx, remove_idx = floes.area[i] < floes.area[j] ?
                        (j, i) : (i, j)
                    fuse_two_floes!(
                        get_floe(floes, keep_idx),
                        get_floe(floes, remove_idx),
                        Δt,
                        floe_settings,
                        prefuse_max_floe_id,
                        rng,
                    )
                    max_floe_id += 1
                    floes.id[keep_idx] = max_floe_id
                end
            end
        end
    end
    return max_floe_id
end

"""
    remove_floes!(
        floes,
        grid,
        domain,
        dissolved,
        floe_settings
    )

Remove floes marked for removal and dissolve floes smaller than minimum floe
area if the dissolve setting is on.
Inputs:
    floes           <StructArray{Floe}> model's floes
    grid            <AbstractRectilinearGrid> model's grid
    domain          <Domain> model's domain
    dissolved       <Matrix{AbstractFloat}> ocean's dissolved field
    floe_settings   <FloeSettings> simulation's settings for making floes
Outputs:
    None. Removes floes that do not continue to the next timestep and reset all
    continuing floes status to active.
"""
function remove_floes!(
    floes,
    grid,
    domain,
    dissolved,
    floe_settings,
)
    for i in reverse(eachindex(floes))
        if floes.status[i].tag != remove && (
            floes.area[i] < floe_settings.min_floe_area ||
            floes.height[i] < floe_settings.min_floe_height
        )
            # Dissolve small/thin floes and add mass to ocean
            dissolve_floe!(
                get_floe(floes, i),
                grid,
                domain,
                dissolved,
            )
            # remove dissolved floes
            StructArrays.foreachfield(
                    field -> deleteat!(field, i),
                    floes,
            )
        elseif floes.status[i].tag == remove
            StructArrays.foreachfield(
                    field -> deleteat!(field, i),
                    floes,
            )
        else  # reset continuing floes' status
            floes.status[i].tag = active
            empty!(floes.status[i].fuse_idx)
        end
    end
    return
end

"""
    simplify_floes!(
        model,
        simp_settings,
        collision_settings,
        floe_settings,
        Δt,
        rng,
    )
Simplify the floe list be smoothing vertices, fusing floes, dissolving floes,
and removing floes as needed. 
Inputs:
    model               <Model> model
    max_floe_id         <Int> maximum floe id in simulation
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    collision_settings  <CollisionSettings> simulation's collision settings
    floe_settings       <FloeSettings>  simulation's settings for making floes
    Δt                  <Int> simulation timestep in seconds
    rng                 <RNG> random number generator
Outputs:
    Updates floe list and removes floe that won't continue to the next timestep
"""
function simplify_floes!(
    model,
    max_floe_id,
    simp_settings,
    collision_settings,
    floe_settings,
    Δt,
    rng,
)
    # Smooth coordinates to reduce total number of vertices
    if simp_settings.smooth_vertices_on
        smooth_floes!(
            model.floes,
            model.domain.topography,
            simp_settings,
            collision_settings,
            floe_settings,
            Δt,
            rng,
        )
    end
    # Fuse floes that have been marked for fusion
    max_floe_id = fuse_floes!(
        model.floes,
        max_floe_id,
        floe_settings,
        Δt,
        rng,
    )

    # Remove floes marked for removal and dissolving
    remove_floes!(
        model.floes,
        model.grid,
        model.domain,
        model.ocean.dissolved,
        floe_settings,
    )
    return max_floe_id
end