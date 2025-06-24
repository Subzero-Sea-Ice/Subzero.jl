"""
Functions needed for calculating bond forces and torques
"""
# Declare initial_bond_endpoints and initial_best_bond_match (minimum bond length is considered when ghosts are involved)
# as a global vector of tuples outside the function
global initial_bond_endpoints = Vector{Tuple{Tuple{Int, Int, Int,Int}, Tuple{Vector{Float64}, Vector{Float64}}}}()
global initial_best_bond_match  = Vector{Tuple{Tuple{Int, Int, Int,Int}, Tuple{Vector{Float64}, Vector{Float64}}}}()

# Function to update hex_bonds when floes are removed from simulation Eg: Open boundariesfunction update_hex_bonds!(floe)
function update_hex_bonds!(floes)
    existing_ids = Set(floes.id)
    #println("existing_ids::", existing_ids )
    for (index, bonds) in enumerate(floes.hex_bonds)
        for j in eachindex(bonds)
            if bonds[j] != 0 && !(bonds[j] in existing_ids)
                floes.hex_bonds[index][j] = 0
            end
        end
    end
end


#Function to calculate differences in velocities at the endpoints of a bond
function get_vel_diff(floes, idx_floe, idx_neighbour, c_floe, c_neighbour)
    u_floe = floes.u[idx_floe] + floes.ξ[idx_floe] * (c_floe[1] - floes.centroid[idx_floe][1])
    v_floe = floes.v[idx_floe] + floes.ξ[idx_floe] * (c_floe[2] - floes.centroid[idx_floe][2])
    u_neighbour = floes.u[idx_neighbour] + floes.ξ[idx_neighbour] * (c_neighbour[1] - floes.centroid[idx_neighbour][1])
    v_neighbour = floes.v[idx_neighbour] + floes.ξ[idx_neighbour] * (c_neighbour[2] - floes.centroid[idx_neighbour][2])
    Δu = u_neighbour - u_floe
    Δv = v_neighbour - v_floe
    return Δu, Δv
end

#Function to calculate midpoint of each face of polygon to find endpoints of a bond
function face_midpoint(p1, p2)
    return [(p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2]
end

#Function to calculate endpoints of a bond
function compute_bond_endpoints(floes, index, idx_floe, bond_id, idx_neighbour_ghost, nfaces)
    bond_idx_floe = findfirst(isequal(bond_id), floes.hex_bonds[index])    
    bond_idx_neighbour = findfirst(isequal(idx_floe), floes.hex_bonds[idx_neighbour_ghost]) 
    c_floe = face_midpoint(
        floes.coords[index][1][bond_idx_floe],
        floes.coords[index][1][mod(bond_idx_floe, nfaces) + 1]
    )
    c_neighbour = face_midpoint(
        floes.coords[idx_neighbour_ghost][1][bond_idx_neighbour],
        floes.coords[idx_neighbour_ghost][1][mod(bond_idx_neighbour, nfaces) + 1]
    )
    return c_floe, c_neighbour
end

# Function to calculate bond length 
function calculate_bond_length(p1, p2)
    return sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
end

# Function to find out the best bond combination between floes (actual+ghost floes) by selecting the combination with minimum distance
function find_closest_best_bond_match(floes, present_bond_endpoints)
    best_bond_match_ordered = Vector{Vector{Tuple{Tuple{Int, Int, Int,Int}, Tuple{Vector{Float64}, Vector{Float64}}}}}()
    seen_combinations = Set{Tuple{Int, Int}}() # Track (idx_floe, bond_id) pairs
    processed_pairs = Set{Tuple{Int, Int}}() # Track processed tuples
    for (index, bonds) in enumerate(floes.hex_bonds)
        idx_floe = floes.id[index]
        for bond_id in bonds
            if bond_id != 0
                key1 = (idx_floe, bond_id)  # The actual order
                key2 = (bond_id, idx_floe)  # The reversed key
                if key1 in processed_pairs || key2 in processed_pairs
                    continue
                end 
                if !(key1 in seen_combinations) && !(key2 in seen_combinations)
                    filtered = filter(x -> (x[1][1], x[1][2]) == key1 || (x[1][1], x[1][2]) == key2, present_bond_endpoints)
                    if !isempty(filtered)
                        min_dist_match = filtered[1]
                        min_dist = calculate_bond_length(filtered[1][2][1], filtered[1][2][2])

                        for match in filtered
                            dist = calculate_bond_length(match[2][1], match[2][2])
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

# Function to find which floes are bonded together (including their ghosts)
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

# Function to extract the bonded floes (To use this in torque calculation)
function extract_best_match(matched_key, key1, key2, found_match)
    if matched_key == key1
        return found_match[2][1], found_match[2][2]
    elseif matched_key == key2
        return found_match[2][2], found_match[2][1]
    end
end

# Function to calculate bond forces and torques at each time step
function timestep_hexagonal_bonds!(
    floes::StructArray{<:Floe{FT}},
    domain,
    consts,
    curr_time,
    hexagonal_bond_settings,
) where {FT<:AbstractFloat}

    nfaces = hexagonal_bond_settings.nfaces
    k_ice  = hexagonal_bond_settings.elastic_constant
    c_ice  = hexagonal_bond_settings.damp_coeff
    global present_bond_endpoints = Vector{Tuple{Tuple{Int, Int, Int, Int}, Tuple{Vector{Float64}, Vector{Float64}}}}()
    global present_best_bond_match = Vector{Tuple{Tuple{Int, Int, Int, Int}, Tuple{Vector{Float64}, Vector{Float64}}}}()
    update_hex_bonds!(floes)   # Updates hex_bonds when floes are removed from simulation
    """
    present_bond_endpoints :: It stores the information of 
    1. Index of floe, 2. Index of neighbour 3. Ghost index of floe 4. Ghost index of neighbour
    4. Endpoints of bonds for all the combinations of bonded floes 
    For Eg: If floe-1 is bonded to floe-4 and both the floes have a ghost too then present_bond_endpoints 
    stores information of endpoints of bonds for all possibilities
    like (floe,neighbour,ghost_floe,ghost_neighbour)=>(1,4,0,0),(1,4,1,0),(1,4,0,1),(1,4,1,1)

    present_best_bond_match :: It stores the information of bonded floes which has least bond length 
    among all combinations stored in present_bond_endpoints. If (1,4,0,1) has least bond length among 
    (1,4,0,0),(1,4,1,0),(1,4,0,1),(1,4,1,1) , then "present_best_bond_match" stores bond endpoints of (1,4,0,1)
    """
    for (index, bonds) in enumerate(floes.hex_bonds)
        idx_floe = floes.id[index]
        idx_floe_ghost = floes.ghost_id[index]
        for bond_id in bonds
            if bond_id != 0
                idx_neighbour_ghost = findall(==(bond_id), floes.id) 
                for idx_neighbour in idx_neighbour_ghost 
                    c_floe, c_neighbour = compute_bond_endpoints(floes,index, idx_floe, bond_id, idx_neighbour, nfaces)
                    push!(present_bond_endpoints, ((idx_floe, bond_id, idx_floe_ghost,floes.ghost_id[idx_neighbour]), (c_floe, c_neighbour)))
                end
            end
        end
    end
    
    if curr_time == 0
        global initial_bond_endpoints = present_bond_endpoints
        global initial_best_bond_match = find_closest_best_bond_match(floes, initial_bond_endpoints)
    end
    present_best_bond_match = find_closest_best_bond_match(floes,present_bond_endpoints) 
    for (index, bonds) in enumerate(floes.hex_bonds)
        fill!(floes.bond_force[index], FT(0))  # Reset both Fx and Fy to zero
        floes.bond_trq[index] = FT(0.0)
        fill!(floes.stress_bond[index], FT(0))  # Reset all stresses to zero
        floe_stress_bond =  zeros(FT, 2, 2)
        # Ensure ind_bond_force is reset at each time step
        if length(floes.ind_bond_force) > index
            for i in eachindex(floes.ind_bond_force[index])
               floes.ind_bond_force[index][i] .= FT(0)  # Set Fx and Fy to zero
            end
        end
        idx_floe = floes.id[index]
        if floes.ghost_id[index] == 0
            for bond_id in bonds
                if bond_id != 0
                    Fx_bond, Fy_bond, torque_xy = 0.0, 0.0, 0.0 # Initialize to avoid scope issues.
                    stress = zeros(FT, 2, 2)
                    key1 = (idx_floe, bond_id) # Actual order
                    key2 = (bond_id, idx_floe) # Reverse order
                    """
                    If  present_best_bond_match stores bond endpoints of (1,4,0,1) {minimum bond length} then
                    find_bond_match() :: present_found_match = stores bond end points (floe and neighbour), present_matched_key = (1,4)
                    ---------------------------------------------------------------------------------------------------------------------
                    extract_best_match() :: Floe is 1, Neighbour is 4 {Index positions of (1,4) matters in (1->floe,4->neighbour,0,1)}
                    c_floe_best_bond_match = stores floe's bond endpoints, c_neighbour_best_bond_match = stores neighbours's bond endpoints based on
                    index positions of floe and neighbour in present_matched_key 
                    """
                    initial_found_match, initial_matched_key = find_bond_match(initial_best_bond_match, key1, key2)
                    present_found_match, present_matched_key = find_bond_match(present_best_bond_match, key1, key2)
                    c_floe_best_bond_match_init, c_neighbour_best_bond_match_init = extract_best_match(initial_matched_key, key1, key2, initial_found_match)
                    c_floe_best_bond_match, c_neighbour_best_bond_match = extract_best_match(present_matched_key, key1, key2, present_found_match)
                    Δx = (c_neighbour_best_bond_match[1] - c_floe_best_bond_match[1]) - (c_neighbour_best_bond_match_init[1] - c_floe_best_bond_match_init[1])
                    Δy = (c_neighbour_best_bond_match[2] - c_floe_best_bond_match[2]) - (c_neighbour_best_bond_match_init[2] - c_floe_best_bond_match_init[2])
                    c_floe = nothing  # Declare and initialize c_floe outside the loop
                    c_neighbour = nothing
                    for (bond_tuple, (found_c_floe, found_c_neighbour)) in present_bond_endpoints
                        if bond_tuple == (idx_floe, bond_id, 0, 0)
                            c_floe = found_c_floe
                            c_neighbour = found_c_neighbour
                            break
                        elseif bond_tuple == (bond_id, idx_floe, 0, 0)
                            c_floe = found_c_neighbour
                            c_neighbour = found_c_floe
                        break
                        end
                    end
                    if c_floe !== nothing && c_neighbour !== nothing
                        index_nei = findfirst(==(bond_id), floes.id)
                        Δu, Δv = get_vel_diff(floes, index, index_nei, c_floe, c_neighbour)
                        # Calculate each bond force for each floe
                        Fx_bond = k_ice * Δx  +  c_ice * Δu
                        Fy_bond = k_ice * Δy  +  c_ice * Δv
                        # Store each Fx_bond,Fy_bond in ind_bond_force
                        idx_nei = findfirst(==(bond_id), bonds)
                        floes.ind_bond_force[index][idx_nei][1] = Fx_bond
                        floes.ind_bond_force[index][idx_nei][2] = Fy_bond
                        # Accumalate all bond forces for each floe
                        floes.bond_force[index][1] += Fx_bond
                        floes.bond_force[index][2] += Fy_bond
                        r_x = c_floe[1] - floes.centroid[index][1]
                        r_y = c_floe[2] - floes.centroid[index][2]
                        # Calculate torque: τ = r × F = r_x * F_y - r_y * F_x for each bond for each floe
                        torque_xy = r_x * Fy_bond - r_y * Fx_bond
                        # Accumulaye all bond torques for each floe
                        floes.bond_trq[index] += torque_xy
                        # Calculate individual bond stresses
                        stress[1,1] = Fx_bond * r_x
                        stress[1,2] = Fx_bond * r_y + Fy_bond * r_x
                        stress[2,2] = Fy_bond * r_y
                        # Accumulate stresses
                        floe_stress_bond[1,1] += stress[1,1]
                        floe_stress_bond[1,2] += stress[1,2]
                        floe_stress_bond[2,2] += stress[2,2]
                      #=   if(curr_time == 20000)
                            step_t = curr_time/20
                            println("Time: $step_t, Floe: $idx_floe, Neighbor: $bond_id, Δx : $Δx, Δy: $Δy")
                            println("Time: $step_t, Floe: $idx_floe, Neighbor: $bond_id, Δu : $Δu, Δv: $Δv")
                            println("Time: $step_t, Floe: $idx_floe, Neighbor: $bond_id, Forces: $Fx_bond, $Fy_bond, Torque: $torque_xy")
                            println("Time: $step_t, Floe: $idx_floe, Neighbor: $bond_id, x: $c_floe")
                            println("Time: $step_t, Floe: $idx_floe, Neighbor: $bond_id", "<cx:>", floes.centroid[index][1], "<cy:>", floes.centroid[index][2])
                            println("Time: $step_t, Floe: $idx_floe, Neighbor: $bond_id, r_x: $r_x, r_y: $r_y")
                            println("Time: $step_t, Floe: $idx_floe, Neighbor: $bond_id, σ_xx: ", stress[1,1], "σ_xy:", stress[1,2],  "σ_yy:", stress[2,2])
                        end  =#
                    else
                        println("c_floe not found for [idx_floe, bond_id, 0, 0]") #add error checking to see if c_floe was not found.
                    end   
                end
            end
        end
        floe_stress_bond[1, 2] *= FT(0.5)
        floe_stress_bond[2, 1] = floe_stress_bond[1, 2]
        floe_stress_bond .*= 1/(floes.area[index] * floes.height[index])
        floes.stress_bond[index] .= floe_stress_bond

        #=  if(curr_time == 20000)
            println("<<Floe:>>",floes.id[index],"<<Force>>", floes.bond_force[index],"<<Torque>>",floes.bond_trq[index])
            println("<<Floe:>>",floes.id[index],"<<Force>>", floes.ind_bond_force[index])
            println("<<Floe:>>",floes.id[index],"<<Bonds>>", floes.hex_bonds[index])
            println("<<Floe:>>",floes.id[index],"<<Area>>", floes.area[index], "<<Height>>", floes.height[index])
            println("<<Floe:>>",floes.id[index],"<<Bond stress>>", floes.stress_bond[index])
        end  =#
    end   
    return
end


