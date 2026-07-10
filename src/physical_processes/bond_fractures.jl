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
    λ_Rb  = bond_settings.λ_Rb

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
                    # 1. Grab local normal and shear forces from history storage
                    Fn_bond = floes.ind_bond_force_nt[idx][edge_idx][1]
                    Ft_bond = floes.ind_bond_force_nt[idx][edge_idx][2]

                    r_i = sqrt(floes.area[idx] / π)
                    r_j = sqrt(floes.area[nei_idx] / π)
                    r_eff = (2 * r_i * r_j) / (r_i + r_j + eps(FT))
                    h_eff = min(floes.height[idx], floes.height[nei_idx])
                    area_b = 2 * λ_Rb * r_eff * h_eff

                    # 3. Calculate scale-dependent critical limits
                    # Normal tensile limit (scaling factor)
                    σt_lim = σt_max_bond * (h_eff / h_mean)^(1.0) * (r_eff / r_mean)^(-0.3)
                    σc_lim = σc_max_bond * (h_eff / h_mean)^(1.0) * (r_eff / r_mean)^(-0.3)
                    τ_lim  = τ_max_bond * (h_eff / h_mean)^(1.0) * (r_eff / r_mean)^(-0.3)
                    F_ten_max   = σt_lim * area_b
                    F_comp_max  = σc_lim * area_b
                    F_shear_max = τ_lim * area_b
                    
                   is_broken = false
            
                   if Fn_bond > 0.0
                       # 1. Tensile Check: Positive Fn means the floes are pulling apart
                       if Fn_bond > F_ten_max
                           is_broken = true
                       end
                   else
                       # 2. Compressive Check: Negative (or zero) Fn means the floes are pushing together
                       if abs(Fn_bond) > F_comp_max
                           is_broken = true
                       end
                   end
                   
                   # 3. Shear Check: Tangential force exceeds the shear capacity
                   if abs(Ft_bond) > F_shear_max
                       is_broken = true
                   end
                    
                    # 5. Break the bond if thresholds are exceeded
                    if is_broken
                        floes.bonds[idx][edge_idx] = 0
                    end
                end
            end
        end
    end

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

