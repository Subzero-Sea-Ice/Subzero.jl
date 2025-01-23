"""
Functions needed for ridging and rafting between floes and boundaries. 
"""

"""
    add_floe_volume!(
        floes,
        idx,
        vol,
        floe_settings,
    )
Add volume to existing floe and update fields that depend on the volume.
Inputs:
    floes           <StructArray{Frloe}> list of floes
    idx             <Int> index of floe to add volume to 
    vol             <AbstractFloat> volume to add to floe
    floe_settings   <FloeSettings> simulation's settings for making floes
Outputs:
    Nothing. Floe's fields are updated to reflect increase in volume.
"""
function add_floe_volume!(
    floes,
    idx,
    vol,
    floe_settings,
)
    # Update floe height, mass, and moment of intertia due to volume change
    init_height = floes.height[idx]
    floes.height[idx] += vol/floes.area[idx]
    if floes.height[idx] > floe_settings.max_floe_height
        floes.height[idx] = floe_settings.max_floe_height
    end
    floes.mass[idx] += vol * floe_settings.ρi
    floes.moment[idx] *= floes.height[idx] / init_height
    # Update all ghost floes of parent
    for gidx in floes.ghosts[idx]
        floes.height[gidx] = floes.height[idx]
        floes.mass[gidx] = floes.mass[idx]
        floes.moment[gidx] = floes.moment[idx]
    end
    return
end

"""
    remove_floe_overlap!(
        floes,
        shrink_idx,
        grow_floe_poly,
        pieces_buffer,
        max_floe_id,
        broken,
        ridgeraft_settings,
        floe_settings,
        rng,  
    )

Removes area/volume of overlap from floe that loses area during ridging/rafting
Inputs:
    floes               <StructArray{Floe}> list of floes
    shrink_idx       <Int> index of floe that loses area
    grow_floe_poly   <Polys> polygon of floe/domain that subsumes area
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    ridgeraft_settings  <RidgeRaftSettings> simulation's ridge/raft settings
    floe_settings       <FloeSettings> simulation's settings for making floes
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    rng                 <AbstractRNG> random number generator
Outputs:
    transfer_vol    <Float> total volume to transfer away from floe
    max_floe_id     <Int> maximum floe id of floe created during overlap removal
    floe_num        <Int> total number of floes created from origianl floe ->
                        one if floe doesn't break, more otherwise
"""
function remove_floe_overlap!(
    floes::StructArray{<:Floe{FT}},
    shrink_idx,
    shrink_parent_idx,
    grow_floe_poly,
    pieces_buffer,
    max_floe_id,
    broken,
    ridgeraft_settings,
    floe_settings,
    simp_settings,
    rng,  
) where {FT <: AbstractFloat}
    # Find new floe shapes and regions
    regions = diff_polys(floes.poly[shrink_idx], grow_floe_poly, FT)
    total_area = zero(FT)
    nregions = 0
    for (i, region) in enumerate(regions)
        regions[i] = simplify_poly(region, simp_settings.tol)
        total_area += GO.area(region)
        nregions += 1
    end
    floe_num = 0  # How many new floes have been created from the regions
    # Changes in area / volume
    transfer_area = floes.area[shrink_idx] - total_area
    transfer_vol = FT(0)
    # If the transfer area is sufficent, create new floes
    if transfer_area > ridgeraft_settings.min_overlap_frac * floes.area[shrink_idx]
        transfer_vol = floes.area[shrink_idx] * floes.height[shrink_idx]
        # Reset shrinking index to parent floe and determine floe shift
        parent_Δx = floes.centroid[shrink_parent_idx][1] -
            floes.centroid[shrink_idx][1]
        parent_Δy = floes.centroid[shrink_parent_idx][2] -
            floes.centroid[shrink_idx][2]
        parent_centroid = floes.centroid[shrink_parent_idx]
        # Update existing floes/ghosts regions
        for region in regions
            region_area = GO.area(region)
            (xmin, xmax), (ymin, ymax) = GI.extent(region)
            Δx, Δy = xmax - xmin, ymax - ymin
            # Region is big enought to be a floe and has okay aspect ratio
            if (
                region_area > floe_settings.min_floe_area &&
                (Δx > Δy ? Δy/Δx : Δx/Δy) > floe_settings.min_aspect_ratio
            )
                floe_num += 1
                new_poly = _translate_poly(FT, region, parent_Δx, parent_Δy)::Polys{FT}
                rmholes!(new_poly)  # remove holes in floe
                new_vol = region_area * floes.height[shrink_idx]
                transfer_vol -= new_vol
                buffer_length = length(pieces_buffer)
                # If this is the first region created, replace original floe
                if floe_num == 1
                    replace_floe!(  # replace parent floe
                        get_floe(floes, shrink_parent_idx),
                        new_poly,
                        new_vol * floe_settings.ρi,
                        floe_settings,
                        rng,
                    )
                    if nregions == 1
                        # if floe doesn't break, propograte changes to ghosts
                        for gidx in floes.ghosts[shrink_parent_idx]
                            # find shift to ghost floe location
                            g_Δx = floes.centroid[gidx][1] - parent_centroid[1]
                            g_Δy = floes.centroid[gidx][2] - parent_centroid[2]
                            # replace ghost floe
                            replace_floe!(
                                get_floe(floes, gidx),
                                new_poly,
                                new_vol * floe_settings.ρi,
                                floe_settings,
                                rng,
                            )
                            # shift ghost floe
                            _translate_floe!(FT, get_floe(floes, gidx), g_Δx, g_Δy)
                        end
                    else
                        # if floe breaks, mark floe and ghosts as broken
                        broken[shrink_parent_idx] = true
                        for gidx in floes.ghosts[shrink_parent_idx]
                            broken[gidx] = true
                            floes.status[gidx].tag = remove
                        end
                        # Update floe identifiers
                        empty!(floes.ghosts[shrink_parent_idx])
                        push!(
                            floes.parent_ids[shrink_parent_idx],
                            floes.id[shrink_parent_idx],
                        )
                        max_floe_id += 1
                        floes.id[shrink_parent_idx] = max_floe_id
                    end
                else  # >1 region, so floe must break and add pieces to buffer
                    push!(
                        pieces_buffer,
                        deepcopy_floe(get_floe(floes, shrink_parent_idx))
                    )
                    buffer_length += 1
                    replace_floe!(
                        get_floe(pieces_buffer, buffer_length),
                        new_poly,
                        new_vol * floe_settings.ρi,
                        floe_settings,
                        rng,
                    )
                    max_floe_id += 1
                    pieces_buffer.id[buffer_length] = max_floe_id
                end
            end
        end
        if floe_num == 0 
            #= If all regions were too small, mark floe for removal and transfer
                all mass to ridging floe
            =#
            floes.status[shrink_parent_idx].tag = remove
        end
    end
    return transfer_vol, max_floe_id, floe_num
end

"""
    floe_floe_ridge!(
        floes,
        idx1,
        idx2,
        floe2,
        overlap_area,
        ridgeraft_settings,
        simp_settings,
        Δt,
        rng
    )
Ridge two floes, updating both in-place and returning any new floes that
resulting from the ridging event.
Inputs:
    floes               <StructArray{Floe}> floe list
    idx1                <Int> index of first floe
    idx2                <Int> index of second floe
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    ridgeraft_settings  <RidgeRaftSettings> simulation's settings for ridging
                            and rafting
    floe_settings       <FloeSettings> simulation's settings for making floes
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    Δt                  <Int> simulation timestep in seconds
    rng                 <RandomNumberGenerator> simulation's random number
                            generator
Outputs:
    Updates floe1 and floe2 and returns any new floes created by ridging
"""
function floe_floe_ridge!(
    floes::StructArray{<:Floe{FT}},
    idx1,
    idx2,
    pieces_buffer,
    max_floe_id,
    broken,
    ridgeraft_settings,
    floe_settings,
    simp_settings,
    Δt,
    rng,
) where {FT}
    # Heights of floes determine which floe subsumes shared area
    f1_h = floes.height[idx1] >= ridgeraft_settings.min_ridge_height
    f2_h = floes.height[idx2] >= ridgeraft_settings.min_ridge_height
    vol = FT(0)
    nregions = 0
    # Determine which floe transfers mass and which gains mass
    gain_mass_idx = 0
    lose_mass_idx = 0 
    if(
        (f1_h && f2_h && rand(rng, FT) >= 1/(1 + (floes.height[idx1]/floes.height[idx2]))) ||
        (f1_h && !f2_h)
    )
        #=
        Either both floes are over min height and we randomly pick floe 1 to
        "subsume" the extra area, or only floe 1 is over min height and it
        gets the extra area
        =#
        gain_mass_idx = idx1
        lose_mass_idx = idx2
    elseif (f1_h && f2_h) ||  (!f1_h && f2_h)
        #=
        Either both floes are over max height and we randomly pick floe 2 to
        "subsumes" the extra area, or only floe 2 is over max height and it
        gets the extra area
        =#
        gain_mass_idx = idx2
        lose_mass_idx = idx1
    end
    if gain_mass_idx > 0 && lose_mass_idx > 0
        # Find parent indices of gain_mass_idx and lose_mass_idx
        lose_parent_idx = lose_mass_idx
        gain_parent_idx = gain_mass_idx
        if floes.ghost_id[lose_mass_idx] != 0
            lose_parent_idx = findfirst(x -> x == floes.id[lose_mass_idx], floes.id)
        end
        if floes.ghost_id[gain_mass_idx] != 0
            gain_parent_idx = findfirst(x -> x == floes.id[gain_mass_idx], floes.id)
        end
        # Inital floe values
        ml, mg = floes.mass[lose_mass_idx], floes.mass[gain_mass_idx]
        Ig = floes.moment[gain_mass_idx]
        xg, yg = floes.centroid[gain_mass_idx]
        # Ridge
        vol, max_floe_id, nregions = remove_floe_overlap!(
            floes,
            lose_mass_idx,
            lose_parent_idx,
            floes.poly[gain_mass_idx],
            pieces_buffer,
            max_floe_id,
            broken,
            ridgeraft_settings,
            floe_settings,
            simp_settings,
            rng,  
        )
        if vol > 0
            add_floe_volume!(
                floes,
                gain_parent_idx,
                vol,
                floe_settings,
            )
            # Conserve momentum
            first_slot = length(pieces_buffer) - nregions + 2
            if nregions < 1
                conserve_momentum_change_floe_shape!(
                    mg, Ig, xg, yg, Δt,
                    get_floe(floes, gain_mass_idx),
                    get_floe(floes, lose_mass_idx),
                )
            elseif nregions == 1
                conserve_momentum_transfer_mass!(floes,
                    lose_mass_idx, gain_mass_idx,
                    ml, mg, Δt,
                )
            else  # floe broke, ghost floes
                conserve_momentum_transfer_mass!(floes,
                    lose_mass_idx, gain_mass_idx,
                    ml, mg, Δt,
                    pieces_buffer, first_slot,
                )
            end
            if !broken[lose_mass_idx]
                update_ghost_timestep_vals!(
                    floes, lose_mass_idx, lose_parent_idx,
                )
            end
            if !broken[gain_mass_idx]
                update_ghost_timestep_vals!(
                    floes, gain_mass_idx, gain_parent_idx,
                )
            end
        end
    end
    return max_floe_id
end

"""
    floe_domain_ridge!(
        floes,
        idx,
        domain_element,
        pieces_buffer,
        max_floe_id,
        broken,
        ridgeraft_settings,
        floe_settings,
        simp_settings,
        Δt,
        rng,
    )

Ridge a floe against a boundary or a topography element and return any excess
floes created by the ridging.
Inputs:
    floes               <StructArray{Floe}> floe list
    idx1                <Int> index of first floe
    domain_element      <AbstractDomainElement> boundary or topography element
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    ridgeraft_settings  <RidgeRaftSettings> simulation's settings for ridge/raft
    floe_settings       <FloeSettings> simulation's settings for making floes
    simp_settings       <SimplificationSettings> simulation's settings for
                            simplification
    Δt                  <Int> simulation timestep in seconds
    rng                 <RandomNumberGenerator> simulation's random number
                            generator
Outputs:
    floe1 is updated with new shape. Return maximum floe id of floes created
"""
function floe_domain_ridge!(
    floes::StructArray{<:Floe{FT}},
    idx,
    domain_element::Union{<:AbstractDomainElement, LazyRow{<:AbstractDomainElement}},
    pieces_buffer,
    max_floe_id,
    broken,
    ridgeraft_settings,
    floe_settings,
    simp_settings,
    Δt,
    rng,
) where {FT}
    # Find parent idx
    parent_idx = idx
    if floes.ghost_id[idx] != 0
        parent_idx = findfirst(x -> x == floes.id[idx], floes.id)
    end
    # Record previous values for momentum conservation
    mass_tmp = floes.mass[idx]
    moment_tmp = floes.moment[idx]
    x_tmp, y_tmp = floes.centroid[idx]
    # Ridge floe with domain element
    vol, max_floe_id, nregions = remove_floe_overlap!(
        floes,
        idx,
        parent_idx,
        domain_element.poly,
        pieces_buffer,
        max_floe_id,
        broken,
        ridgeraft_settings,
        floe_settings,
        simp_settings,
        rng,  
    )
    if vol > 0 && nregions > 0
        # Determine if extra volume should be added to floe or to the domain
        if rand(rng, FT) > ridgeraft_settings.domain_gain_probability
            current_slot = length(pieces_buffer) - nregions + 2
            tot_area = floes.area[idx] + sum(pieces_buffer.area[current_slot:end])
            for i in 1:nregions
                if i == 1
                    region_frac = floes.area[idx] / tot_area
                    add_floe_volume!(
                        floes,
                        idx,
                        vol * region_frac,
                        floe_settings,
                    )
                else
                    region_frac = pieces_buffer.area[current_slot] / tot_area
                    add_floe_volume!(
                        pieces_buffer,
                        current_slot,
                        vol * region_frac,
                        floe_settings,
                    )
                    current_slot += 1
                end
            end
        end
        # Update floe velocities to conserve momentum as domain element has no
        # momentum, but floe now has less mass if ridge was successful
        if nregions == 1
            conserve_momentum_change_floe_shape!(
                mass_tmp,
                moment_tmp,
                x_tmp,
                y_tmp,
                Δt,
                get_floe(floes, idx),
            )
        end
        if !broken[idx]
            update_ghost_timestep_vals!(floes, idx, parent_idx)
        end
    end
    return max_floe_id
end

"""
    floe_floe_raft!(
        floes,
        idx1,
        idx2,
        pieces_buffer,
        max_floe_id,
        broken,
        ridgeraft_settings,
        floe_settings,
        simp_settings,
        Δt,
        rng,
    )

Raft two floes, updating both in-place and returning any new floes that
resulting from the rafting event.
Inputs:
    floes               <StructArray{Floe}> floe list
    idx1                <Int> index of first floe
    idx2                <Int> index of second floe
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    ridgeraft_settings  <RidgeRaftSettings> simulation's ridge/raft settings
    floe_settings       <FloeSettings> simultion's settings for making floes
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    Δt                  <Int> simulation timestep in seconds
    rng                 <RandomNumberGenerator> simulation's random number
                            generator
Outputs:
    Updates floe1 and floe2 and returns any new floes created by rafting
"""
function floe_floe_raft!(
    floes::StructArray{<:Floe{FT}},
    idx1,
    idx2,
    pieces_buffer,
    max_floe_id,
    broken,
    ridgeraft_settings,
    floe_settings,
    simp_settings,
    Δt,
    rng,
) where {FT}
    vol = FT(0)
    nregions = 0
    # Based on height ratio, pick which floe subsumes shares area
    # Default is floe 2 subsumes mass from floe 1
    gain_mass_idx = idx2
    lose_mass_idx = idx1
    if rand(rng, FT) >= 1/(1 + (floes.height[idx1]/floes.height[idx2]))
        # Floe 1 subsumes mass from floe 2
        gain_mass_idx = idx1
        lose_mass_idx = idx2
    end
    # Find parents
    lose_parent_idx = lose_mass_idx
    gain_parent_idx = gain_mass_idx
    if floes.ghost_id[lose_parent_idx] != 0
        lose_parent_idx = findfirst(x -> x == floes.id[lose_mass_idx], floes.id)
    end
    if floes.ghost_id[gain_mass_idx] != 0
        gain_parent_idx = findfirst(x -> x == floes.id[gain_mass_idx], floes.id)
    end
    # Inital floe values
    ml, mg = floes.mass[lose_mass_idx], floes.mass[gain_mass_idx]
    Ig = floes.moment[gain_mass_idx]
    xg, yg = floes.centroid[gain_mass_idx]
    # Raft
    vol, max_floe_id, nregions = remove_floe_overlap!(
        floes,
        lose_mass_idx,
        lose_parent_idx,
        floes.poly[gain_mass_idx],
        pieces_buffer,
        max_floe_id,
        broken,
        ridgeraft_settings,
        floe_settings,
        simp_settings,
        rng,  
    )
    if vol > 0 && nregions > 0
        # Add extra area/volume to floe 2
        add_floe_volume!(
            floes,
            gain_parent_idx,
            vol,
            floe_settings,
        )
        first_slot = length(pieces_buffer) - nregions + 2
        if nregions == 0
            conserve_momentum_change_floe_shape!(
                mg, Ig, xg, yg, Δt,
                get_floe(floes, gain_mass_idx),
                get_floe(floes, lose_mass_idx),
            )
        elseif nregions == 1
            conserve_momentum_transfer_mass!(floes,
                lose_mass_idx, gain_mass_idx,
                ml, mg, Δt,
            )
        else  # floe broke, ghost floes
            conserve_momentum_transfer_mass!(floes,
                lose_mass_idx, gain_mass_idx,
                ml, mg, Δt,
                pieces_buffer, first_slot,
            )
        end
        if !broken[lose_mass_idx]
            update_ghost_timestep_vals!(floes, lose_mass_idx, lose_parent_idx)
        end
        if !broken[gain_mass_idx]
            update_ghost_timestep_vals!(floes, gain_mass_idx, gain_parent_idx)
        end
    end
    return max_floe_id
end

"""
    floe_domain_raft!(
        floes,
        idx1,
        domain_element,
        pieces_buffer,
        max_floe_id,
        broken,
        ridgeraft_settings,
        floe_settings,
        simp_settings,
        Δt
        rng,
    )

Raft a floe against a boundary or a topography element and return any excess
floes created by the rafting. This is equivalent to ridging.
Inputs:
    floes               <StructArray{Floe}> floe list
    idx1                <Int> index of first floe
    domain_element      <AbstractDomainElement> boundary or topography element
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    ridgeraft_settings  <RidgeRaftSettings> ridge/raft settings
    floe_settings       <FloeSettings> simulation's settings for making floes
    simp_settings       <SimplificationSettings> simplification settings
    Δt                  <Int> simulation timestep in seconds
    rng                 <RandomNumberGenerator> simulation's random number
                            generator
Outputs:
    floe1 is updated with new shape. If any new floes are created by rafting
    they are returned, else nothing.
"""
floe_domain_raft!(
    floes,
    idx1,
    domain_element::Union{AbstractDomainElement, LazyRow{<:AbstractDomainElement}},
    pieces_buffer,
    max_floe_id,
    broken,
    ridgeraft_settings,
    floe_settings,
    simp_settings,
    Δt,
    rng,
) = floe_domain_ridge!(
    floes,
    idx1,
    domain_element,
    pieces_buffer,
    max_floe_id,
    broken,
    ridgeraft_settings,
    floe_settings,
    simp_settings,
    Δt,
    rng,
)

"""
    timestep_ridging_rafting!(
        floes,
        pieces_buffer,
        domain,
        max_floe_id,
        ridgeraft_settings::RidgeRaftSettings{FT},
        floe_settings
        simp_settings,
        Δt,
        rng,
    )

Ridge and raft floes that meet probability and height criteria.
Inputs:
    floes               <StructArray{Floe}> simulation's list of floes
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    domain              <Domain> simulation's domain
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    ridgeraft_settings  <RidgeRaftSettings> ridge/raft settings
    floe_settings       <FloeSettings> simulation's settings for making floes
    simp_settings       <SimplificationSettings> simplification settings
    Δt                  <Int> length of timestep in seconds
    rng                 <RandomNumberGenerator> simulation's rng
Outputs:
    Updates floes post ridging and rafting and adds any new pieces to the pieces
    buffer to be made into new floes.
"""
function timestep_ridging_rafting!(
    floes,
    pieces_buffer,
    domain,
    max_floe_id,
    ridgeraft_settings::RidgeRaftSettings{FT},
    floe_settings,
    simp_settings,
    Δt,
    rng = Xoshiro(),
) where {FT <: AbstractFloat}
    broken = fill(false, length(floes))
    interactions_list = zeros(Int, size(floes.interactions[1], 1))
    for i in eachindex(floes)
        #=
            Floe is active in sim, hasn't broken, interacted with other floes,
            isn't too thick, and meets probability check to ridge or raft
        =#
        ridge = floes.height[i] <= ridgeraft_settings.max_floe_ridge_height &&
            rand(rng, FT) <= ridgeraft_settings.ridge_probability
        raft = floes.height[i] <= ridgeraft_settings.max_floe_raft_height &&
            rand(rng, FT) <= ridgeraft_settings.raft_probability
        if (
            (ridge || raft) &&
            floes.status[i].tag == active &&
            !broken[i] &&
            floes.num_inters[i] > 0
        )
            ninters = 0
            # Find unique interactions that haven't already been calculated
            for row in 1:floes.num_inters[i]
                j = Int(floes.interactions[i][row, floeidx])
                min_area = j > 0 ?
                    min(floes.area[i], floes.area[j]) :  # floe-floe interaction
                    floes.area[i] # floe-domain interaction

                # floes/domain overlap (not ghost interaction copied to parent)
                valid_interaction = false
                if i < j && !broken[j] && floes.status[j].tag == active
                    valid_interaction |= potential_interaction(
                        floes.centroid[i], floes.centroid[j],
                        floes.rmax[i], floes.rmax[j],
                    )
                elseif j == -1
                    valid_interaction |=
                        abs(floes.centroid[i][2] - domain.north.val) <
                        floes.rmax[i]
                elseif j == -2
                    valid_interaction |=
                        abs(floes.centroid[i][2] - domain.south.val) <
                        floes.rmax[i]
                elseif j == -3
                    valid_interaction |=
                        abs(floes.centroid[i][1] - domain.east.val) <
                        floes.rmax[i]
                elseif j == -4
                    valid_interaction |=
                        abs(floes.centroid[i][1] - domain.west.val) <
                        floes.rmax[i]
                elseif j < 0
                    valid_interaction |= potential_interaction(
                        floes.centroid[i],
                        domain.topography.centroid[-(j + 4)],
                        floes.rmax[i],
                        domain.topography.rmax[-(j + 4)],
                    )
                end
                if (valid_interaction && !(j in interactions_list) &&
                    1e-6 < floes.interactions[i][row, overlap]/min_area < 0.95
                )
                    ninters += 1
                    if ninters <= length(interactions_list)
                        interactions_list[ninters] = j
                    else
                        push!(interactions_list, j)
                    end
                end     
            end 
            for k in 1:ninters
                # if floe-floe interaction
                j = interactions_list[k]
                if j > 0 && !broken[i] && !broken[j]
                    if (  # if both floes aren't too thick -> ridge
                        ridge &&
                        floes.height[i] <= ridgeraft_settings.max_floe_ridge_height &&
                        floes.height[j] <= ridgeraft_settings.max_floe_ridge_height
                    )
                        max_floe_id = floe_floe_ridge!(
                            floes,
                            i,
                            j,
                            pieces_buffer,
                            max_floe_id,
                            broken,
                            ridgeraft_settings,
                            floe_settings,
                            simp_settings,
                            Δt,
                            rng,
                        )
                    elseif ( # if both floes aren't too thick -> raft
                        raft &&
                        floes.height[i] <= ridgeraft_settings.max_floe_raft_height &&
                        floes.height[j] <= ridgeraft_settings.max_floe_raft_height
                    )
                        max_floe_id = floe_floe_raft!(
                            floes,
                            i,
                            j,
                            pieces_buffer,
                            max_floe_id,
                            broken,
                            ridgeraft_settings,
                            floe_settings,
                            simp_settings,
                            Δt,
                            rng,
                        )
                    end
                elseif j < 0 && !broken[i]  # if floe-domain interaction
                    if (  # if floe isn't too thick -> ridge
                        ridge &&
                        floes.height[i] <= ridgeraft_settings.max_domain_ridge_height
                    )
                        max_floe_id = floe_domain_ridge!(
                            floes,
                            i,
                            get_domain_element(domain, j),
                            pieces_buffer,
                            max_floe_id,
                            broken,
                            ridgeraft_settings,
                            floe_settings,
                            simp_settings,
                            Δt,
                            rng,
                        )

                    elseif (  # if floe isn't too thick -> raft
                        raft &&
                        floes.height[i] <= ridgeraft_settings.max_domain_raft_height
                    )
                        max_floe_id = floe_domain_raft!(
                            floes,
                            i,
                            get_domain_element(domain, j),
                            pieces_buffer,
                            max_floe_id,
                            broken,
                            ridgeraft_settings,
                            floe_settings,
                            simp_settings,
                            Δt,
                            rng,
                        )
                    end
                end
            end
        end
    end
    return max_floe_id
end
