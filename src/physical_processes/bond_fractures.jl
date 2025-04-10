"""
Function to implement bond breakage.
"""
function bond_fracture_floes!(
    floes::StructArray{<:Floe{FT}},
    hexagonal_bond_fracture_settings,
) where {FT<:AbstractFloat}
    
    σc_bond = hexagonal_bond_fracture_settings.σc_bond
    rad_bond = hexagonal_bond_fracture_settings.rad_bond
    Ft = σc_bond * (π * (rad_bond^2))  # Critical force threshold
    #println("<<Bonds in fracture>>",floes.hex_bonds)
    for (index, bonds) in enumerate(floes.hex_bonds)
        #println("<<Bonds Force>> Floe:",index, "<<-->>",floes.ind_bond_force[index])
        idx_floe = floes.id[index]
        for bond_id in bonds
            if bond_id != 0
                idx_nei = findfirst(==(bond_id), bonds)
                idx_floe_in_nei = findfirst(==(idx_floe), floes.hex_bonds[bond_id])
                Fx_bond, Fy_bond = floes.ind_bond_force[index][idx_nei]
                resul_force = sqrt(Fx_bond^2 + Fy_bond^2)
                if resul_force >= Ft
                    floes.hex_bonds[index][idx_nei] = 0            # Remove bond in floe
                    floes.hex_bonds[bond_id][idx_floe_in_nei] = 0  # Remove bond in neighbour
                end
            end
        end
    end
end



