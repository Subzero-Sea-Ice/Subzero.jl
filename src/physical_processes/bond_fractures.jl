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
    #println("<<Bonds Force in time step>>",length(floes.hex_bonds))
    for (index, bonds) in enumerate(floes.hex_bonds)
        #println("<<Bonds Force>> Floe:",index, "<<-->>",floes.ind_bond_force[index])
        idx_floe = floes.id[index]
        for bond_id in bonds
            if bond_id != 0
                idx_neighbour = findfirst(==(bond_id), floes.id)
                idx_nei = findfirst(==(bond_id), bonds)
                idx_floe_in_nei = findfirst(==(idx_floe), floes.hex_bonds[idx_neighbour])
                Fx_bond, Fy_bond = floes.ind_bond_force[index][idx_nei]
                resul_force = sqrt(Fx_bond^2 + Fy_bond^2)
                if resul_force >= Ft
                    floes.hex_bonds[index][idx_nei] = 0            # Remove bond in floe
                    floes.hex_bonds[idx_neighbour][idx_floe_in_nei] = 0  # Remove bond in neighbour
                   #=  # Make sure that the ghost floe's hex_bonds are also set to zero
                    floe_has_ghost = any((floes.id[i] == idx_floe && floes.ghost_id[i] == 1) for i in eachindex(floes.id))
                    floe_nei_has_ghost = any((floes.id[i] == bond_id && floes.ghost_id[i] == 1) for i in eachindex(floes.id))

                    frac_floe_ghost_idx = floe_has_ghost ? findall(==(idx_floe), floes.id) : Int[]
                    frac_floe_nei_ghost_idx = floe_nei_has_ghost ? findall(==(bond_id), floes.id) : Int[]

                    if !isempty(frac_floe_ghost_idx)
                        for floe_ghost_idx in frac_floe_ghost_idx
                            floes.hex_bonds[floe_ghost_idx][idx_nei] = 0
                        end
                    end
                    if !isempty(frac_floe_nei_ghost_idx)
                        for floe_nei_ghost_idx in frac_floe_nei_ghost_idx
                            floes.hex_bonds[floe_nei_ghost_idx][idx_floe_in_nei] = 0
                        end
                    end =#

                end
            end
        end
    end
end



