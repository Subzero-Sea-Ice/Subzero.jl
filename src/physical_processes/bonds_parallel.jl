"""
Functions needed for calculating bond forces and torques
"""
function update_bonds!(floes)
    existing_ids = Set(floes.id)
    #println("existing_ids::", existing_ids )
    for (index, bonds) in enumerate(floes.bonds)
        for j in eachindex(bonds)
            if bonds[j] != 0 && !(bonds[j] in existing_ids)
                floes.bonds[index][j] = 0
            end
        end
    end
end

function polygon_edges(poly_coords::Vector{Vector{Float64}})
    n = length(poly_coords)
    pts = [GB.Point2(c[1], c[2]) for c in poly_coords]
    [GB.Line(pts[i], pts[mod1(i+1, n)]) for i in 1:n]
end

function midpoint(edge::GB.Line)
    p1, p2 = edge.points
    GB.Point2((p1[1] + p2[1])/2, (p1[2] + p2[2])/2)
end

# Function to check if point is inside polygon (outer ring)
function point_in_polygon(p::GB.Point2, poly_coords::Vector{Vector{Float64}})
    x, y = p[1], p[2]
    n = length(poly_coords)
    inside = false
    j = n
    tol = 1e-9

    for i in 1:n
        xi, yi = poly_coords[i]
        xj, yj = poly_coords[j]

        # --- Check if point lies on this edge ---
        cross = (x - xi) * (yj - yi) - (y - yi) * (xj - xi)
        if abs(cross) < tol
            if (min(xi, xj) - tol ≤ x ≤ max(xi, xj) + tol) &&
               (min(yi, yj) - tol ≤ y ≤ max(yi, yj) + tol)
                return true  # point is on the boundary
            end
        end

        # --- Ray casting for inside/outside test ---
        intersect = ((yi > y) != (yj > y)) &&
                    (x < (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)

        inside ⊻= intersect
        j = i
    end

    return inside
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
function timestep_bonds!(
    floes::StructArray{<:Floe{FT}},
    domain,
    consts,
    curr_time,
    bond_settings,
) where {FT<:AbstractFloat}

    λ_Rb  = bond_settings.λ_Rb
    λ_lb  = bond_settings.λ_lb
    ζ_b   = bond_settings.ζ_b
    SF_b  = bond_settings.SF_b
    k_b  = bond_settings.k_b
    c_b  = bond_settings.c_b

    # --- Step 1. Update bonds due to fracturing ---
    update_bonds!(floes)

    # --- Step 2. Reset per-edge forces (parallel safe) ---
    Threads.@threads for id in eachindex(floes)
        for i in eachindex(floes.ind_bond_force[id])
            floes.ind_bond_force[id][i] .= FT(0)
        end
    end

    # --- Step 3. Build ID lookup (serial, small overhead) ---
    id_to_indices = Dict{Int, Vector{Int}}()
    for (i, f) in enumerate(floes)
        push!(get!(id_to_indices, f.id, Int[]), i)
    end

    # --- Step 4. Main parallel bond computation ---
    Threads.@threads for idx in eachindex(floes)
        tid = Threads.threadid()
        id = floes.id[idx]
        if floes.ghost_id[idx] == 0  # process only parent floes
            floes.bond_force[idx] .= FT(0)
            edges_i = polygon_edges(floes.coords[idx][1])
            num_edges_i = length(edges_i) - 1

            for edge_idx in 1:num_edges_i
                neighbor = floes.bonds[idx][edge_idx]

                # only compute once (id < neighbor)
                if neighbor != 0 && id < neighbor
                    # Get all ghost indices (differs geometrically)
                    i_indices = id_to_indices[id]
                    neighbor_indices = get(id_to_indices, neighbor, Int[])

                    best_dist = Inf
                    best_dx = 0.0
                    best_dy = 0.0
                    best_pair = nothing

                    # Search across all (i_idx, neighbor_idx) combinations
                    for i_idx in i_indices
                        edges_i2 = polygon_edges(floes.coords[i_idx][1])
                        edge_i = edges_i2[edge_idx]
                        mid_i = midpoint(edge_i)

                        for neighbor_idx in neighbor_indices
                            edges_n = polygon_edges(floes.coords[neighbor_idx][1])
                            nei_edge_idx = findfirst(x -> x == id, floes.bonds[neighbor_idx])
                            edge_n = edges_n[nei_edge_idx]
                            mid_n = midpoint(edge_n)

                            dx = mid_n[1] - mid_i[1]
                            dy = mid_n[2] - mid_i[2]
                            dist = sqrt(dx^2 + dy^2)

                            #if point_in_polygon(mid_i, floes.coords[neighbor_idx][1]) ||
                            #   point_in_polygon(mid_n, floes.coords[i_idx][1])
                            #    dx = dy = dist = 0.0
                            #end

                            if dist < best_dist
                                best_dist = dist
                                best_dx, best_dy = dx, dy
                                best_pair = (i_idx, edge_idx, neighbor_idx, nei_edge_idx, mid_i, mid_n)
                            end
                        end
                    end

                    # --- Assign forces for best pair ---
                    if best_pair !== nothing
                        i_idx2, edge_idx2, neighbor_idx, nei_edge_idx, mid_i2, mid_n2 = best_pair
                        thread_id = Threads.threadid()
                        parent_id = floes.id[i_idx2]
                        floe_idx_act = findfirst(x -> x == parent_id, floes.id) 
                        nei_idx_act = findfirst(x -> x == neighbor, floes.id)
    
                        Δu, Δv = get_vel_diff(floes, i_idx2, neighbor_idx, mid_i2, mid_n2)
                        r_i = sqrt(floes.area[i_idx2] / π)
                        r_j = sqrt(floes.area[neighbor_idx] / π)
                        r_eff = (2 * r_i * r_j) / (r_i + r_j)
                        h_eff = min(floes.height[i_idx2], floes.height[neighbor_idx])
                        area_b = 2 * λ_Rb * r_eff * h_eff

                        dxc = floes.centroid[neighbor_idx][1] - floes.centroid[i_idx2][1]
                        dyc = floes.centroid[neighbor_idx][2] - floes.centroid[i_idx2][2]
                        len_b = sqrt(dxc^2 + dyc^2)
                        k_b = consts.E * area_b / len_b
                        m_b = (floes.mass[i_idx2] * floes.mass[neighbor_idx]) /
                              (floes.mass[i_idx2] + floes.mass[neighbor_idx])
                        c_b = 2 * ζ_b * sqrt(k_b * m_b)       

                        Fx_bond = best_dx * k_b + c_b * Δu
                        Fy_bond = best_dy * k_b + c_b * Δv

                        floes.ind_bond_force[floe_idx_act][edge_idx2][1] = Fx_bond
                        floes.ind_bond_force[floe_idx_act][edge_idx2][2] = Fy_bond
    
                        #val_x = floes.ind_bond_force[floe_idx_act][edge_idx2][1]
                        #val_y = floes.ind_bond_force[floe_idx_act][edge_idx2][2]
			#if curr_time == 100
			#    println("Thread $thread_id | Bond $(parent_id) -> $(nei_idx_act) | Target Index $(floe_idx_act) (Edge $(edge_idx2))")
			#    println(" $(parent_id) -> $(nei_idx_act) [Geometric Inputs]: best_dx = $best_dx, best_dy = $best_dy (best_dist = $best_dist)")
			#    println(" $(parent_id) -> $(nei_idx_act) [Kinematic Inputs]: Δu = $Δu, Δv = $Δv")
			#    println(" $(parent_id) -> $(nei_idx_act) [Property Inputs]: k_b = $k_b, c_b = $c_b (mass_eff = $m_b)")
			#    println(" $(parent_id) -> $(nei_idx_act) [CALC]: Fx = $Fx_bond, Fy = $Fy_bond")
			#    println(" $(parent_id) -> $(nei_idx_act) [WRITTEN]: X = $val_x, Y = $val_y")
			#    println(" $(parent_id) -> $(nei_idx_act) [INDEX CHECK]: Winner Floe Index: $i_idx2 vs. Target Index: $floe_idx_act (Must be same if geometry matters!)")
			#    println(" $(parent_id) -> $(nei_idx_act) [INDEX CHECK]: Winner Nei Index: $neighbor_idx vs. Target Index: $nei_idx_act (Must be same if geometry matters!)")
			#end

                    end
                end
            end
        end
    end

    # --- Step 4a. Mirror forces for id > neighbor ---
    for idx in eachindex(floes)
        id = floes.id[idx]
        bonds_i = floes.bonds[idx]
        num_edges_i = length(bonds_i)
        if floes.ghost_id[idx] == 0
            for edge_idx in 1:num_edges_i
               neighbor = bonds_i[edge_idx]
               if neighbor != 0 && id > neighbor
                     nei_idx_act = findfirst(x -> x == neighbor, floes.id)
                     nei_edge_idx = findfirst(x -> x == id, floes.bonds[nei_idx_act])
                     floes.ind_bond_force[idx][edge_idx][1] = -floes.ind_bond_force[nei_idx_act][nei_edge_idx][1]
                     floes.ind_bond_force[idx][edge_idx][2] = -floes.ind_bond_force[nei_idx_act][nei_edge_idx][2]
               end
            end
        end
    end

    # --- Step 4b. Compute total bond force for every floe (after all threads finish) ---
    Threads.@threads for idx in eachindex(floes)
        sx = zero(FT)
        sy = zero(FT)
        for f in floes.ind_bond_force[idx]
            sx += f[1]
            sy += f[2]
        end
        floes.bond_force[idx][1] = sx
        floes.bond_force[idx][2] = sy
    end

    # --- Step 5. Parallel torque & stress computation ---
    Threads.@threads for idx in eachindex(floes)
        floes.bond_trq[idx] = FT(0)
        fill!(floes.stress_bond[idx], FT(0))
        floe_stress_bond = zeros(FT, 2, 2)

        if floes.ghost_id[idx] == 0
            edges_i = polygon_edges(floes.coords[idx][1])
            num_edges_i = length(edges_i) - 1
            for edge_idx in 1:num_edges_i
                edge_i = edges_i[edge_idx]
                mid_i = midpoint(edge_i)
                r_x = mid_i[1] - floes.centroid[idx][1]
                r_y = mid_i[2] - floes.centroid[idx][2]
                Fx_bond, Fy_bond = floes.ind_bond_force[idx][edge_idx]
                torque_xy = r_x * Fy_bond - r_y * Fx_bond

                floes.bond_trq[idx] += torque_xy
                floe_stress_bond[1,1] += Fx_bond * r_x
                floe_stress_bond[1,2] += Fx_bond * r_y + Fy_bond * r_x
                floe_stress_bond[2,2] += Fy_bond * r_y
            end
        end

        floe_stress_bond[1,2] *= FT(0.5)
        floe_stress_bond[2,1] = floe_stress_bond[1,2]
        denom = floes.area[idx] * floes.height[idx] + eps(FT)
        floe_stress_bond .*= 1 / denom
        floes.stress_bond[idx] .= floe_stress_bond
    end

    return
end

