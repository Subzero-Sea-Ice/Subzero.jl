"""
Function to implement bond breakage.
"""
function bond_fracture_floes!(
    floes::StructArray{<:Floe{FT}},
    bond_fracture_settings,
    bond_settings,
    ) where {FT<:AbstractFloat}
    σt_max_bond = bond_fracture_settings.σt_max_bond
    σc_max_bond = bond_fracture_settings.σc_max_bond
    τ_max_bond  = bond_fracture_settings.τ_max_bond

  
    #println("<<Bonds Force in time step>>",length(floes.bonds))
    for (index, bonds) in enumerate(floes.bonds)
        #println("<<Bonds Force>> Floe:",index, "<<-->>",floes.ind_bond_force[index])
        idx_floe = floes.id[index]
        for bond_id in bonds
            if bond_id != 0
                idx_neighbour = findfirst(==(bond_id), floes.id)
                idx_nei = findfirst(==(bond_id), bonds)
                idx_floe_in_nei = findfirst(==(idx_floe), floes.bonds[idx_neighbour])
                Fx_bond, Fy_bond = floes.ind_bond_force[index][idx_nei]
                resul_force = sqrt(Fx_bond^2 + Fy_bond^2)
                r_i = sqrt(floes.area[index] / π)
                r_j = sqrt(floes.area[idx_neighbour] / π)
                r_eff = (2 * r_i * r_j) / (r_i + r_j)                          # reduced (harmonic) radius
                h_eff = min(floes.height[index], floes.height[idx_neighbour])   # conservative effective height
                area_b = 2 * bond_settings.λ_Rb * r_eff * h_eff                # bond cross-section area
                Ft = σt_max_bond * area_b
                if resul_force >= Ft
                    floes.bonds[index][idx_nei] = 0                  # Remove bond in floe
                    floes.bonds[idx_neighbour][idx_floe_in_nei] = 0  # Remove bond in neighbour   
                end
            end
        end
    end
end
