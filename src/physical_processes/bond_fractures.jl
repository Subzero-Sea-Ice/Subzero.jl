"""
Function to implement bond breakage.
"""
function polygon_edges(poly_coords::Vector{Vector{Float64}})
    n = length(poly_coords)
    pts = [GB.Point2(c[1], c[2]) for c in poly_coords]
    [GB.Line(pts[i], pts[mod1(i+1, n)]) for i in 1:n]
end

function bond_fracture_floes!(
    floes::StructArray{<:Floe{FT}},
    bond_fracture_settings,
    bond_settings,
) where {FT<:AbstractFloat}

    σt_max_bond = bond_fracture_settings.σt_max_bond
    σc_max_bond = bond_fracture_settings.σc_max_bond
    τ_max_bond  = bond_fracture_settings.τ_max_bond
    c0          = bond_fracture_settings.cohesion
    μ_bond      = bond_fracture_settings.mu_bond
    h_mean      = bond_fracture_settings.h_mean
    r_mean      = bond_fracture_settings.r_mean
    λ_Rb        = bond_settings.λ_Rb

    # --- Step 1. Build fast ID lookup table ---
    id_to_indices = Dict{Int, Int}()
    for (i, f) in enumerate(floes)
        id_to_indices[f.id] = i
    end

    # --- Step 2. Parallel fracture checks ---
    Threads.@threads for idx in eachindex(floes)
        id = floes.id[idx]
        
        if floes.ghost_id[idx] == 0
            bonds_list = floes.bonds[idx]
            edges_i = polygon_edges(floes.coords[idx][1]) 
            num_edges = length(bonds_list)
            
            for edge_idx in 1:num_edges
                neighbor_id = bonds_list[edge_idx]
                
                # Symmetrical safety check: process each bond exactly once
                if neighbor_id != 0 && id < neighbor_id
                    
                    nei_idx = id_to_indices[neighbor_id]
                    edges_j = polygon_edges(floes.coords[nei_idx][1])
                    floei_edge = edges_i[edge_idx]
                    edge_idx_j = findfirst(x -> x == id, floes.bonds[nei_idx])
                    floej_edge = edges_j[edge_idx_j]
                    
                    # Midline averaging calculations to find physical edge length (Lb)
                    p1 = (floei_edge.points[1] + floej_edge.points[2]) / 2.0  
                    p2 = (floei_edge.points[2] + floej_edge.points[1]) / 2.0  
                    Lb = norm(p2 - p1)

                    # 1. Grab local normal force, shear force, and bending moment from history storage
                    Fn_bond = floes.ind_bond_force_nt[idx][edge_idx][1]
                    Ft_bond = floes.ind_bond_force_nt[idx][edge_idx][2]
                    Mb_bond = floes.ind_bond_bm[idx][edge_idx]

                    # 2. Geometric calculations
                    r_i = sqrt(floes.area[idx] / π)
                    r_j = sqrt(floes.area[nei_idx] / π)
                    r_eff = (2 * r_i * r_j) / (r_i + r_j + eps(FT))
                    h_eff = min(floes.height[idx], floes.height[nei_idx])
                    area_b = Lb * h_eff

                    # 3. Calculate scale-dependent critical stress limits
                    σt_lim = σt_max_bond * (h_eff / h_mean)^(0.5) * (r_eff / r_mean)^(-0.5)
                    σc_lim = σc_max_bond * (h_eff / h_mean)^(0.5) * (r_eff / r_mean)^(-0.5)
                    τ_lim  = τ_max_bond  * (h_eff / h_mean)^(0.5) * (r_eff / r_mean)^(-0.5)

                    # 4. Stress Calculations
                    # Normal axial stress (tension is positive, compression is negative)
                    σ_axial = Fn_bond / (area_b + eps(FT))

                    # Bending stress
                    σ_bend = (6.0 * abs(Mb_bond)) / (h_eff * Lb^2 + eps(FT))

                    # Interface shear stress
                    τ_shear = abs(Ft_bond) / (area_b + eps(FT))
       
                    # Check maximum possible tension and compression for combined axial and bending
                    peak_tension = σ_axial +  σ_bend           
                    peak_compression = -σ_axial +  σ_bend
                    break_mode = 0

                    # Gate 1 :  Tensile/ Shear Ellipse Failure
                    if peak_tension > 0.0
                        failure_index = (peak_tension/σt_lim)^2 + (τ_shear/τ_lim)^2
                        if failure_index >= 1.0
                            break_mode = 1
                        end
                    end

                    # Gate 2 :  Compressive Crushing
                    if break_mode == 0 && peak_compression > 0.0
                        if peak_compression >= σc_lim
                            break_mode = 2
                        end
                    end

                    # 5. Break the bond if any thresholds are exceeded
                    if break_mode > 0
                        edge_idx_j = findfirst(x -> x == id, floes.bonds[nei_idx])
                        floes.bond_break_mode[idx][edge_idx] = break_mode
                        floes.bond_break_mode[nei_idx][edge_idx_j] = break_mode
                        floes.bonds[idx][edge_idx] = 0
                    end
                end
            end
        end
    end

    # --- Step 3. Symmetrical bond clean-up check ---
    for (index, bonds) in enumerate(floes.bonds)
        floe_id = floes.id[index]
        for bond_id in bonds
            if bond_id != 0
                 idx_neighbour = findfirst(==(bond_id), floes.id)
                 idx_pos = findfirst(==(bond_id), bonds)
                 bonds_nei = floes.bonds[idx_neighbour] 
                 if findfirst(==(floe_id), bonds_nei) === nothing
                    floes.bonds[index][idx_pos] = 0
                 end 
            end
        end
    end

end
