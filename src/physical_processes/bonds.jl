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
    update_bonds!(floes)     # To take exiting floes into account
   

    # --- Step 2. Reset per-edge global xy forces (Keep nt historical states intact!) ---
    Threads.@threads for id in eachindex(floes)
        for i in eachindex(floes.ind_bond_force_xy[id])
            floes.ind_bond_force_xy[id][i] .= FT(0)
            floes.ind_bond_stiff_nt[id][i] .= FT(0)
        end
    end

    # --- Step 3. Build ID lookup table for al floes, including ghosts---
   id_to_indices = Dict{Int, Vector{Int}}()
   for (i, f) in enumerate(floes)
       push!(get!(id_to_indices, f.id, Int[]), i)
   end
    # println("<<No. of floes>>",length(floes.id))
    # println("<<ids>>",floes.id)
    # println("<<Bonds>>",floes.bonds)
    # println("<<Breakmode>>",floes.bond_break_mode)
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
                    # Find all ghost floes of current floe (geometry differs, that's why!!)
                    i_indices = id_to_indices[id]
                    neighbor_indices = get(id_to_indices, neighbor_id, Int[])
                    best_dist = Inf
                    best_dx = 0.0
                    best_dy = 0.0
                    best_pair = nothing
                    for i_idx in i_indices # Go through all floes (parent + ghost)
                        edges_i2 = polygon_edges(floes.coords[i_idx][1])
                        edge_i = edges_i2[edge_idx]
                        mid_i = midpoint(edge_i)
                        for neighbor_idx in neighbor_indices  # Go through all neighbouring floes (parent + ghost)
                            edges_n = polygon_edges(floes.coords[neighbor_idx][1])
                            nei_edge_idx = findfirst(x -> x == id, floes.bonds[neighbor_idx])
                            edge_n = edges_n[nei_edge_idx]
                            mid_n = midpoint(edge_n)
                            dx = mid_n[1] - mid_i[1]
                            dy = mid_n[2] - mid_i[2]
                            dist = sqrt(dx^2 + dy^2)
                            if dist < best_dist
                                best_dist = dist
                                best_dx, best_dy = dx, dy
                                best_pair = (i_idx, edge_idx, neighbor_idx, nei_edge_idx, mid_i, mid_n)
                            end
                            # best_pair = (index pos of floe, index pos of edge of floe connected to neighbour floe,
                            #             index pos of neighbour floe, edge of nei floe connected to floe )
                            # 1 -> Floe10 -> [20,3,4,5] # 2 -> Floe20 -> [5,6,7,10] best_pair:(1,1,2,4)
                        end
                    end

                    if best_pair !== nothing
                        idx2, edge_idx2, nei_idx2, nei_edge_idx2, mid_i2, mid_n2 = best_pair
                        # Midline averaging calculations
                        # id -> idx(act of parent) -> idx2 (parent/ghost)
                        # neighbor_id -> nei_idx2
                        floei_edge = edges_i[edge_idx2]
                        edges_j = polygon_edges(floes.coords[nei_idx2][1])
                        floej_edge = edges_j[nei_edge_idx2]
                        p1 = (floei_edge.points[1] + floej_edge.points[2]) / 2.0  
                        p2 = (floei_edge.points[2] + floej_edge.points[1]) / 2.0  
                        m = (p1 + p2) / 2.0
                        Lb = norm(p2 - p1)
      
                        # Thread-safe local variables for unit vectors
                        local_tx, local_ty = (p2[1] - p1[1]) / (Lb + eps(FT)), (p2[2] - p1[2]) / (Lb + eps(FT))
                        local_nx, local_ny = -local_ty, local_tx
                    
                        # Outward pointing alignment check via centroid vector 
                        ci = floes.centroid[idx2]
                        cj = floes.centroid[nei_idx2]
                        R_ij = cj - ci
                        if (local_nx * R_ij[1] + local_ny * R_ij[2]) < 0
                            local_nx = -local_nx
                            local_ny = -local_ny
                        end
                    
                        v_rel_x, v_rel_y = get_vel_diff(floes, idx2, nei_idx2, m)
                        v_rel_n  = v_rel_x * local_nx + v_rel_y * local_ny
                        v_rel_t  = v_rel_x * local_tx + v_rel_y * local_ty
                        # Geometrical properties
                        r_i = sqrt(floes.area[idx2] / π)
                        r_j = sqrt(floes.area[nei_idx2] / π)
                        r_eff = (2 * r_i * r_j) / (r_i + r_j + eps(FT))
                        h_eff = min(floes.height[idx2], floes.height[nei_idx2])
                        area_b = Lb * h_eff
                        dxc = cj[1] - ci[1]
                        dyc = cj[2] - ci[2]
                        len_b = sqrt(dxc^2 + dyc^2) + eps(FT)
                        I_z = (h_eff * Lb^3) / 12.0
                    
                        # Stiffness criteria
                        k_b_val = consts.E * area_b / (len_b + eps(FT))
                        G_b_val = (consts.E / (2 * (1 + consts.ν))) * area_b / (len_b + eps(FT))
                        k_theta = consts.E * I_z / (len_b + eps(FT))

                        # 1. Load previous values from historical storage arrays
                        Fn_old = floes.ind_bond_force_nt[idx][edge_idx][1]
                        Ft_old = floes.ind_bond_force_nt[idx][edge_idx][2]
                        nx_old = floes.ind_bond_normal_old[idx][edge_idx][1]
                        ny_old = floes.ind_bond_normal_old[idx][edge_idx][2]
                        if nx_old == 0.0 && ny_old == 0.0
                            nx_old = local_nx
                            ny_old = local_ny
                        end
                        cos_dθ = nx_old * local_nx + ny_old * local_ny
                        sin_dθ = nx_old * local_ny - ny_old * local_nx
                        Fn_old_rotated =  Fn_old * cos_dθ + Ft_old * sin_dθ
                        Ft_old_rotated = -Fn_old * sin_dθ + Ft_old * cos_dθ 
                        m_b = (floes.mass[idx2] * floes.mass[nei_idx2]) /
                                  (floes.mass[idx2] + floes.mass[nei_idx2] + eps(FT))
                        c_b = 2 * ζ_b * sqrt(k_b_val * m_b)
                        Fn_bond = (k_b_val * v_rel_n * Δt) + (Fn_old_rotated * ζ_b)
                        Ft_bond = (G_b_val * v_rel_t * Δt) + (Ft_old_rotated * ζ_b)
                        floes.ind_bond_normal_old[idx][edge_idx][1] = local_nx
                        floes.ind_bond_normal_old[idx][edge_idx][2] = local_ny
                        # if id < neighbor_id  
                        #     println(
                        #         "Time: $curr_time | Bond [$id <-> $neighbor_id] | Edge: $edge_idx\n" *
                        #         "  -> Normal Vec: Old=[$(round(nx_old, digits=4)), $(round(ny_old, digits=4))] " *
                        #         "--> New=[$(round(local_nx, digits=4)), $(round(local_ny, digits=4))]\n" *
                        #         "  -> Rotation  : cos(dθ)=$(round(cos_dθ, digits=4)), sin(dθ)=$(round(sin_dθ, digits=4))\n" *
                        #         "  -> Forces    : Fn_old=$(round(Fn_old, digits=2)) -> Fn_rot=$(round(Fn_old_rotated, digits=2)) -> Fn_new=$(round(Fn_bond, digits=2))"
                        #     )
                        # end
                        # Convert accumulated local states back to global x-y frames
                        Fx_bond = (Fn_bond * local_nx) + (Ft_bond * local_tx)
                        Fy_bond = (Fn_bond * local_ny) + (Ft_bond * local_ty)
                    
                       # Store current accumulated states back onto Floe-i
                        floes.ind_bond_force_nt[idx][edge_idx][1] = Fn_bond
                        floes.ind_bond_force_nt[idx][edge_idx][2] = Ft_bond   
                        floes.ind_bond_force_xy[idx][edge_idx][1] = Fx_bond
                        floes.ind_bond_force_xy[idx][edge_idx][2] = Fy_bond 
                        
                        # Cache active material stiffnesses for diagnostics
                        floes.ind_bond_stiff_nt[idx][edge_idx][1] = k_b_val
                        floes.ind_bond_stiff_nt[idx][edge_idx][2] = G_b_val
   
                        # Store bond areas
                        floes.ind_bond_prop[idx][edge_idx][1] = area_b
                        floes.ind_bond_prop[idx][edge_idx][2] = I_z
                        floes.ind_bond_prop[idx][edge_idx][3] = Lb
                        # Calculating bending moment
                        Δω_z = floes.ξ[nei_idx2] - floes.ξ[idx2]
                        M_old = floes.ind_bond_bm[idx][edge_idx]
                        M_bond = (k_theta * Δω_z * Δt) + (M_old * ζ_b)
                        floes.ind_bond_bm[idx][edge_idx] = M_bond

                       #if curr_time == 0
                       #        println("==================================================")
                       #        println("Bond propertiesat  curr_time = 0")
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
                     floes.ind_bond_bm[idx][edge_idx] = -floes.ind_bond_bm[nei_idx_act][nei_edge_idx]
                     floes.ind_bond_prop[idx][edge_idx][1] = floes.ind_bond_prop[nei_idx_act][nei_edge_idx][1]
                     floes.ind_bond_prop[idx][edge_idx][2] = floes.ind_bond_prop[nei_idx_act][nei_edge_idx][2]
                     floes.ind_bond_prop[idx][edge_idx][3] = floes.ind_bond_prop[nei_idx_act][nei_edge_idx][3]
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
                BM_bond = floes.ind_bond_bm[idx][edge_idx]
                floes.bond_trq[idx] += (torque_xy + BM_bond)
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
