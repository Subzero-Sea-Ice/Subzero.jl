"""
Functions needed for calculating bond forces and torques
"""

#function calc_bond_force(
#    
#) where {FT<:AbstractFloat}
#return
#end

#function calc_torque(
#    
#    ) where {FT<:AbstractFloat}
#return
#end

function midpoint(p1, p2)
    return [(p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2]
end

function timestep_hexagonal_bonds!(
    floes::StructArray{<:Floe{FT}},
    domain,
    consts,
    curr_time,
    hexagonal_bond_settings,
) where {FT<:AbstractFloat}
   if curr_time == 0
    coord_int = floes.coords    
   end
   
   for (index, polygon) in enumerate(floes.coords)
    bonds = floes.hex_bonds[index]
    idx_floe = floes.id[index]  
    for bond_id in bonds          
        if bond_id != 0
            idx_neighbour = bond_id          
            bond_idx_floe = findfirst(isequal(bond_id), bonds)          
            bonds_neighbour = floes.hex_bonds[idx_neighbour]
            bond_idx_neighbour = findfirst(isequal(idx_floe), bonds_neighbour) 
            p1_floe,p2_floe = floes.coords[idx_floe][1][bond_idx_floe],floes.coords[idx_floe][1][mod(bond_idx_floe, 6) + 1] 
            c_floe = midpoint(p1_floe,p2_floe)   
            p1_neighbour,p2_neighbour = floes.coords[idx_neighbour][1][bond_idx_neighbour],floes.coords[idx_neighbour][1][mod(bond_idx_neighbour, 6) + 1]     
            c_neighbour = midpoint(p1_neighbour,p2_neighbour)    
        end
    end
   end
end
