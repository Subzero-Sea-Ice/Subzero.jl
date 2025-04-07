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

    for (index, bonds) in enumerate(floes.hex_bonds)
        for bond_id in bonds
            if bond_id != 0
                idx_nei = findfirst(==(bond_id), bonds)
                Fx_bond, Fy_bond = floes.ind_bond_force[index][idx_nei]
                resul_force = sqrt(Fx_bond^2 + Fy_bond^2)
                if resul_force >= Ft
                    floes.hex_bonds[index][idx_nei] = 0  # Remove bond
                end
            end
        end
    end
end



