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
function get_vel_diff(floes, idx_floe, idx_neighbour, mid)
    u_floe = floes.u[idx_floe] - floes.ξ[idx_floe] * (mid[2] - floes.centroid[idx_floe][2])
    v_floe = floes.v[idx_floe] + floes.ξ[idx_floe] * (mid[1] - floes.centroid[idx_floe][1])
    u_neighbour = floes.u[idx_neighbour] - floes.ξ[idx_neighbour] * (mid[2] - floes.centroid[idx_neighbour][2])
    v_neighbour = floes.v[idx_neighbour] + floes.ξ[idx_neighbour] * (mid[1] - floes.centroid[idx_neighbour][1])
    Δu = u_neighbour - u_floe
    Δv = v_neighbour - v_floe
    return Δu, Δv
end

function timestep_bonds!(
    floes::StructArray{<:Floe{FT}},
    domain,
    consts,
    Δt,
    curr_time,
    bond_settings,
) where {FT<:AbstractFloat}

    λ_Rb  = bond_settings.λ_Rb
    λ_lb  = bond_settings.λ_lb
    ζ_b   = bond_settings.ζ_b
    SF_b  = bond_settings.SF_b
    
    # --- Step 1. Update bonds due to fracturing ---
    update_bonds!(floes)

    # --- Step 2. Reset per-edge global xy forces (Keep nt historical states intact!) ---
    Threads.@threads for id in eachindex(floes)
        for i in eachindex(floes.ind_bond_force_xy[id])
            floes.ind_bond_force_xy[id][i] .= FT(0)
            floes.ind_bond_stiff_nt[id][i] .= FT(0)
        end
    end

    # --- Step 3. Build ID lookup table ---
    id_to_indices = Dict{Int, Int}()
    for (i, f) in enumerate(floes)
        id_to_indices[f.id] = i
    end

    # --- Step 5. Main parallel bond computation ---
    Threads.@threads for idx in eachindex(floes)
        id = floes.id[idx]
        
        if floes.ghost_id[idx] == 0  # Process parent elements
            floes.bond_force[idx] .= FT(0)
            edges_i = polygon_edges(floes.coords[idx][1]) 
            num_edges_i = length(edges_i) - 1
            
            for edge_idx in 1:num_edges_i
                neighbor_id = floes.bonds[idx][edge_idx]
                
                # Symmetrical safety: process bond interaction exactly once
                if neighbor_id != 0 && id < neighbor_id
                    
                    nei_idx = id_to_indices[neighbor_id]
                    edges_j = polygon_edges(floes.coords[nei_idx][1])
                    
                    ci = floes.centroid[idx]
                    cj = floes.centroid[nei_idx]
                    
                    floei_edge = edges_i[edge_idx]
                    
                    # Find matching shared edge index position on our neighbor
                    edge_idx_j = findfirst(x -> x == id, floes.bonds[nei_idx])
                    floej_edge = edges_j[edge_idx_j]
                    
                    # Midline averaging calculations
                    p1 = (floei_edge.points[1] + floej_edge.points[2]) / 2.0  
                    p2 = (floei_edge.points[2] + floej_edge.points[1]) / 2.0  
                    m = (p1 + p2) / 2.0
                    Lb = norm(p2 - p1)
                    
                    # Thread-safe local variables for unit vectors
                    local_tx, local_ty = (p2[1] - p1[1]) / (Lb + eps(FT)), (p2[2] - p1[2]) / (Lb + eps(FT))
                    local_nx, local_ny = -local_ty, local_tx
                    
                    # Outward pointing alignment check via Centroid vector compass
                    R_ij = cj - ci
                    if (local_nx * R_ij[1] + local_ny * R_ij[2]) < 0
                        local_nx = -local_nx
                        local_ny = -local_ny
                    end
                    
                    v_rel_x, v_rel_y = get_vel_diff(floes, idx, nei_idx, m)
                    v_rel_n  = v_rel_x * local_nx + v_rel_y * local_ny
                    v_rel_t  = v_rel_x * local_tx + v_rel_y * local_ty
                    
                    # Geometrical properties
                    r_i = sqrt(floes.area[idx] / π)
                    r_j = sqrt(floes.area[nei_idx] / π)
                    r_eff = (2 * r_i * r_j) / (r_i + r_j + eps(FT))
                    h_eff = min(floes.height[idx], floes.height[nei_idx])
                    area_b = 2 * λ_Rb * r_eff * h_eff
                    dxc = floes.centroid[nei_idx][1] - floes.centroid[idx][1]
                    dyc = floes.centroid[nei_idx][2] - floes.centroid[idx][2]
                    len_b = sqrt(dxc^2 + dyc^2) + eps(FT)
                    
                    # Stiffness criteria
                    k_b_val = consts.E * area_b / (len_b + eps(FT))
                    G_b_val = (consts.E / (2 * (1 + consts.ν))) * area_b / (len_b + eps(FT))

                    # 1. Load previous values from historical storage arrays
                    Fn_old = floes.ind_bond_force_nt[idx][edge_idx][1]
                    Ft_old = floes.ind_bond_force_nt[idx][edge_idx][2]
                    omega_bond = 0.5 * (floes.ξ[idx] + floes.ξ[nei_idx])
                    delta_theta = omega_bond * Δt
                    cos_dθ = cos(delta_theta)
                    sin_dθ = sin(delta_theta)
                    Fn_old_rotated =  Fn_old * cos_dθ + Ft_old * sin_dθ
                    Ft_old_rotated = -Fn_old * sin_dθ + Ft_old * cos_dθ 
                    Fn_bond = (k_b_val * v_rel_n * Δt) + (Fn_old_rotated * ζ_b)
                    Ft_bond = (G_b_val * v_rel_t * Δt) + (Ft_old_rotated * ζ_b)
                    #println("Time: $curr_time | Bond [$id <-> $neighbor_id] | Old: ($(round(Fn_old, digits=2)), $(round(Ft_old, digits=2))) -> Rotated: ($(round(Fn_old_rotated, digits=2)), $(round(Ft_old_rotated, digits=2))) | New Incr: ($(round(k_b_val * v_rel_n * Δt, digits=2)), $(round(G_b_val * v_rel_t * Δt, digits=2)))")
                    #println("Time: $curr_time | dθ (rad): $delta_theta | Old: $(Fn_old) -> Rotated: $(Fn_old_rotated)")
                    # Convert accumulated local states back to global x-y frames
                    Fx_bond = (Fn_bond * local_nx) + (Ft_bond * local_tx)
                    Fy_bond = (Fn_bond * local_ny) + (Ft_bond * local_ty)
                    
                    # Store current accumulated states back onto Floe I
                    floes.ind_bond_force_nt[idx][edge_idx][1] = Fn_bond
                    floes.ind_bond_force_nt[idx][edge_idx][2] = Ft_bond   
                    floes.ind_bond_force_xy[idx][edge_idx][1] = Fx_bond
                    floes.ind_bond_force_xy[idx][edge_idx][2] = Fy_bond 
                    
                    # Cache active material stiffnesses for diagnostics
                    floes.ind_bond_stiff_nt[idx][edge_idx][1] = k_b_val
                    floes.ind_bond_stiff_nt[idx][edge_idx][2] = G_b_val
                    #if curr_time == 50
                    #        println("==================================================")
                    #        println("BOND DIAGNOSTICS AT curr_time = 0")
                    #        println("Interaction Pair: Floe ID ($id) <-> Neighbor ID ($neighbor_id)")
                    #        println("--------------------------------------------------")
                    #        println("  p1         : ", p1)
                    #        println("  p2         : ", p2)
                    #        println("  Midpoint   : ", m)
                    #        println("  Length (Lb): ", Lb)
                    #        println("  Normal Vec : nx = ", local_nx, ", ny = ", local_ny)
                    #        println("  Tangential : tx = ", local_tx, ", ty = ", local_ty)
                    #end
                    #println("Time: $curr_time | All local bond forces: ", floes.ind_bond_force_nt)
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
                     floes.ind_bond_force_xy[idx][edge_idx][1] = -floes.ind_bond_force_xy[nei_idx_act][nei_edge_idx][1]
                     floes.ind_bond_force_xy[idx][edge_idx][2] = -floes.ind_bond_force_xy[nei_idx_act][nei_edge_idx][2]
                     floes.ind_bond_force_nt[idx][edge_idx][1] = -floes.ind_bond_force_nt[nei_idx_act][nei_edge_idx][1]
                     floes.ind_bond_force_nt[idx][edge_idx][2] = -floes.ind_bond_force_nt[nei_idx_act][nei_edge_idx][2]
                     floes.ind_bond_stiff_nt[idx][edge_idx][1] = floes.ind_bond_stiff_nt[nei_idx_act][nei_edge_idx][1]
                     floes.ind_bond_stiff_nt[idx][edge_idx][2] = floes.ind_bond_stiff_nt[nei_idx_act][nei_edge_idx][2]
               end
            end
        end
    end 

    # --- Step 4b. Compute total bond force for every floe (after all threads finish) ---
    Threads.@threads for idx in eachindex(floes)
        sx = zero(FT)
        sy = zero(FT)
        for f in floes.ind_bond_force_xy[idx]
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
                Fx_bond, Fy_bond = floes.ind_bond_force_xy[idx][edge_idx]
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
