"""
Function to implement bond breakage
"""
function bond_fracture_floes!(
    floes::StructArray{<:Floe{FT}},
    hexagonal_bond_settings,
) where {FT<:AbstractFloat}
σc_bond = hexagonal_bond_settings.σc_bond
rad_bond  = hexagonal_bond_settings.rad_bond
Ft = σc_bond*(π*(rad_bond^2))
#println(Ft)
end


