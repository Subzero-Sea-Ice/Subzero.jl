"""
Functions needed for calculating bond forces and torques
"""
# Declare initial_midpoints as a global vector of tuples outside the function
global initial_bond_endpoints = Vector{Tuple{Tuple{Int, Int, Int,Int}, Tuple{Vector{Float64}, Vector{Float64}}}}()
global initial_best_bond_match  = Vector{Tuple{Tuple{Int, Int, Int,Int}, Tuple{Vector{Float64}, Vector{Float64}}}}()

function midpoint(p1, p2)
    return [(p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2]
end


function compute_bond_midpoints(floes, index, idx_floe, idx_neighbour, idx_neighbour_ghost, nfaces)
    bond_idx_floe = findfirst(isequal(idx_neighbour), floes.hex_bonds[idx_floe])
    bond_idx_neighbour = findfirst(isequal(idx_floe), floes.hex_bonds[idx_neighbour])

    c_floe = midpoint(
        floes.coords[index][1][bond_idx_floe],
        floes.coords[index][1][mod(bond_idx_floe, nfaces) + 1]
    )

    c_neighbour = midpoint(
        floes.coords[idx_neighbour_ghost][1][bond_idx_neighbour],
        floes.coords[idx_neighbour_ghost][1][mod(bond_idx_neighbour, nfaces) + 1]
    )

    return c_floe, c_neighbour
end

function calculate_distance(p1, p2)
    return sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
end

function find_closest_best_bond_match(floes, present_bond_endpoints)
    best_bond_match_ordered = Vector{Vector{Tuple{Tuple{Int, Int, Int,Int}, Tuple{Vector{Float64}, Vector{Float64}}}}}()
    seen_combinations = Set{Tuple{Int, Int}}() # Track (idx_floe, bond_id) pairs
    processed_pairs = Set{Tuple{Int, Int}}() # Track processed tuples
    for (index, bonds) in enumerate(floes.hex_bonds)
        idx_floe = floes.id[index]
        for bond_id in bonds
            if bond_id != 0
                key1 = (idx_floe, bond_id)
                key2 = (bond_id, idx_floe)  # The reversed key
                if key1 in processed_pairs || key2 in processed_pairs
                    continue
                end 
                if !(key1 in seen_combinations) && !(key2 in seen_combinations)
                    filtered = filter(x -> (x[1][1], x[1][2]) == key1 || (x[1][1], x[1][2]) == key2, present_bond_endpoints)
                    if !isempty(filtered)
                        min_dist_match = filtered[1]
                        min_dist = calculate_distance(filtered[1][2][1], filtered[1][2][2])

                        for match in filtered
                            dist = calculate_distance(match[2][1], match[2][2])
                            if dist < min_dist
                                min_dist = dist
                                min_dist_match = match
                            end
                        end
                        push!(best_bond_match_ordered, [min_dist_match])
                        push!(seen_combinations, key1)
                        push!(seen_combinations, key2)
                    end
                end
                push!(processed_pairs, key1)
                push!(processed_pairs, key2)
            end
        end
    end
    return best_bond_match_ordered
end

function find_bond_match(best_bond_match, key1, key2)
    found_match = nothing
    matched_key = nothing
    for match_vector in best_bond_match
        if !isempty(match_vector)
            match_key = (match_vector[1][1][1], match_vector[1][1][2])
            # Check both orderings and store the matched key
            if match_key == key1
                found_match = match_vector[1]
                matched_key = key1  # Store which key matched
                break
            elseif match_key == key2
                found_match = match_vector[1]
                matched_key = key2  # Store which key matched
                break
            end
        end
    end
    return found_match, matched_key
end

function extract_best_match(matched_key, key1, key2, found_match)
    if matched_key == key1
        return found_match[2][1], found_match[2][2]
    elseif matched_key == key2
        return found_match[2][2], found_match[2][1]
    end
end


function timestep_hexagonal_bonds!(
    floes::StructArray{<:Floe{FT}},
    domain,
    consts,
    curr_time,
    hexagonal_bond_settings,
) where {FT<:AbstractFloat}

    nfaces = hexagonal_bond_settings.nfaces
    k_ice  = hexagonal_bond_settings.elastic_constant

    global present_bond_endpoints = Vector{Tuple{Tuple{Int, Int, Int, Int}, Tuple{Vector{Float64}, Vector{Float64}}}}()
    global present_best_bond_match = Vector{Tuple{Tuple{Int, Int, Int, Int}, Tuple{Vector{Float64}, Vector{Float64}}}}()

    for (index, bonds) in enumerate(floes.hex_bonds)
        idx_floe = floes.id[index]
        idx_floe_ghost = floes.ghost_id[index]
        for bond_id in bonds
            if bond_id != 0
                idx_neighbour_ghost = findall(==(bond_id), floes.id)   
                for idx_neighbour in idx_neighbour_ghost 
                    c_floe, c_neighbour = compute_bond_midpoints(floes,index, idx_floe, bond_id, idx_neighbour, nfaces)
                    push!(present_bond_endpoints, ((idx_floe, bond_id, idx_floe_ghost,floes.ghost_id[idx_neighbour]), (c_floe, c_neighbour)))
                end
            end
        end
    end
    if curr_time==0
        global initial_bond_endpoints = present_bond_endpoints
        global initial_best_bond_match = find_closest_best_bond_match(floes, initial_bond_endpoints)
    end
    present_best_bond_match = find_closest_best_bond_match(floes,present_bond_endpoints)
 
    for (index, bonds) in enumerate(floes.hex_bonds)
        fill!(floes.bond_force[index], FT(0))  # Reset both Fx and Fy to zero
        floes.bond_trq[index] = FT(0.0)
        idx_floe = floes.id[index]
        idx_floe_ghost = floes.ghost_id[index]
        if floes.ghost_id[index] == 0
            for bond_id in bonds
                if bond_id != 0
                    Fx_bond, Fy_bond, torque_xy = 0.0, 0.0, 0.0 # Initialize to avoid scope issues.
                    key1 = (idx_floe, bond_id) # Actual order
                    key2 = (bond_id, idx_floe) # Reverse order
                    initial_found_match, initial_matched_key = find_bond_match(initial_best_bond_match, key1, key2)
                    present_found_match, present_matched_key = find_bond_match(present_best_bond_match, key1, key2)
                    c_floe_best_bond_match_init, c_neighbour_best_bond_match_init = extract_best_match(initial_matched_key, key1, key2, initial_found_match)
                    c_floe_best_bond_match, c_neighbour_best_bond_match = extract_best_match(present_matched_key, key1, key2, present_found_match)
                    Δx = (c_neighbour_best_bond_match[1] - c_floe_best_bond_match[1]) - (c_neighbour_best_bond_match_init[1] - c_floe_best_bond_match_init[1])
                    Δy = (c_neighbour_best_bond_match[2] - c_floe_best_bond_match[2]) - (c_neighbour_best_bond_match_init[2] - c_floe_best_bond_match_init[2])
                    k_ice = 1.0 # Example value, adjust accordingly
                    Fx_bond = k_ice * Δx
                    Fy_bond = k_ice * Δy
                    floes.bond_force[index][1] += Fx_bond
                    floes.bond_force[index][2] += Fy_bond
                    c_floe = nothing  # Declare and initialize c_floe outside the loop
                    for (bond_tuple, (found_c_floe, c_neighbour)) in present_bond_endpoints
                        if bond_tuple == (idx_floe, bond_id, 0, 0)
                            c_floe = found_c_floe
                            break
                        elseif bond_tuple == (bond_id, idx_floe, 0, 0)
                            c_floe = c_neighbour
                        break
                        end
                    end
                    if c_floe !== nothing #check if c_floe was found.
                        r_x = c_floe[1] - floes.centroid[index][1]
                        r_y = c_floe[2] - floes.centroid[index][1]
                        torque_xy = r_x * Fy_bond - r_y * Fx_bond
                        #println("Time: $curr_time, Floe: $idx_floe, Neighbor: $bond_id, Forces: $Fx_bond, $Fy_bond, Torque: $torque_xy")
                        #println("Time: $curr_time, Floe: $idx_floe, Neighbor: $bond_id, r_x: $r_x, r_y: $r_y")
                        # Calculate torque: τ = r × F = r_x * F_y - r_y * F_x
                        floes.bond_trq[index] += torque_xy
                    else
                        println("c_floe not found for [idx_floe, bond_id, 0, 0]") #add error checking to see if c_floe was not found.
                    end   
                end
            end
        end
        #println(floes.id[index],"<<Force>>", floes.bond_force[index],"<<Torque>>",floes.bond_trq[index])
    end   

    return
end
