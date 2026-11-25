"""
Function to implement bond breakage.

Bond encoding in `floes.bonds[idx][k]`:
    > 0  id of the bonded neighbour (active bond)
      0  never bonded / no neighbour
     -1  broken in tension
     -2  broken in compression
     -3  broken in shear
"""

# Reset one edge after its bond has broken: the bond entry becomes -mode and the per-edge
# history is cleared, so a broken bond cannot leave stale forces, moments or normals behind.
@inline function _clear_bond!(floes, idx, k, mode)
    floes.bonds[idx][k] = -mode
    fill!(floes.ind_bond_force_nt[idx][k],   0)
    fill!(floes.ind_bond_force_xy[idx][k],   0)
    fill!(floes.ind_bond_stiff_nt[idx][k],   0)
    fill!(floes.ind_bond_normal_old[idx][k], 0)
    fill!(floes.ind_bond_prop[idx][k],       0)
    floes.ind_bond_bm[idx][k] = 0
    return
end

# True if stress x has reached the limit `lim`, within a relative tolerance that absorbs round-off.
# Only positive stresses can break a bond.
@inline _reaches(x, lim, tol) = x > 0 && x >= lim * (1 - tol)

function bond_fracture_floes!(
    floes::StructArray{<:Floe{FT}},
    bond_fracture_settings,
    bond_settings;          # unused, kept so existing callers do not change
    rtol::Real = 100 * eps(FT),     # relative tolerance on the stress limits (round-off)
) where {FT<:AbstractFloat}

    σt_max = FT(bond_fracture_settings.σt_max_bond)
    σc_max = FT(bond_fracture_settings.σc_max_bond)
    τ_max  = FT(bond_fracture_settings.τ_max_bond)
    inv_h_mean = 1 / FT(bond_fracture_settings.h_mean)
    inv_r_mean = 1 / FT(bond_fracture_settings.r_mean)
    ϵ = eps(FT)
    tol = FT(rtol)

    # --- Step 1. Lookup tables, built once ---
    id_to_indices = Dict{Int, Vector{Int}}()    # id -> parent + ghosts
    parent_idx    = Dict{Int, Int}()            # id -> parent only
    for i in eachindex(floes)
        push!(get!(id_to_indices, floes.id[i], Int[]), i)
        floes.ghost_id[i] == 0 && (parent_idx[floes.id[i]] = i)
    end

    # --- Step 2. Parallel fracture checks ---
    # Each bond is checked once, from the floe with the smaller id. This step only WRITES to the
    # checked floe's own slot of `pend`, so it is thread safe. Bonds are not changed yet because
    # step 3 still needs the neighbour id to find the partner edge.
    pend = Vector{Vector{NTuple{2, Int}}}(undef, length(floes))   # (edge, mode) per floe, lazily created
    Threads.@threads for idx in eachindex(floes)
        floes.ghost_id[idx] == 0 || continue
        id      = floes.id[idx]
        bonds_i = floes.bonds[idx]
        for k in eachindex(bonds_i)
            nid = bonds_i[k]
            (nid <= 0 || nid < id) && continue          # only active bonds, each checked once
            j = get(parent_idx, nid, 0)
            j == 0 && continue

            # forces and bond properties from history storage
            fnt    = floes.ind_bond_force_nt[idx][k]
            Fn     = fnt[1]
            Ft     = fnt[2]
            Mb     = floes.ind_bond_bm[idx][k]
            prop   = floes.ind_bond_prop[idx][k]
            area_b = prop[1]
            Lb     = prop[3]

            # scale-dependent critical stress limits (same factor for all three limits)
            r_i   = sqrt(floes.area[idx] / FT(π))
            r_j   = sqrt(floes.area[j]   / FT(π))
            r_eff = (2 * r_i * r_j) / (r_i + r_j + ϵ)
            h_eff = min(floes.height[idx], floes.height[j])
            scale = (h_eff * inv_h_mean) / sqrt(r_eff * inv_r_mean)   # (h/h̄)^1 * (r/r̄)^(-1/2)
            σt_lim = σt_max * scale
            σc_lim = σc_max * scale
            τ_lim  = τ_max  * scale

            # stresses: axial (tension > 0), bending, interface shear
            σ_axial = Fn / (area_b + ϵ)
            σ_bend  = (6 * abs(Mb)) / (h_eff * Lb^2 + ϵ)
            τ_shear = abs(Ft) / (area_b + ϵ)
            peak_tension     =  σ_axial + σ_bend
            peak_compression = -σ_axial + σ_bend

            mode = 0
            if _reaches(peak_tension, σt_lim, tol)                      # 1: tensile
                mode = 1
            elseif _reaches(τ_shear, τ_lim, tol)                        # 3: shear
                mode = 3
            elseif _reaches(peak_compression, σc_lim, tol)              # 2: compressive crushing
                mode = 2
            end

            if mode > 0
                isassigned(pend, idx) || (pend[idx] = NTuple{2, Int}[])
                push!(pend[idx], (k, mode))
            end
        end
    end

    # --- Step 3. Apply breaks on BOTH floes of each bond (serial, only floes that had a break) ---
    changed = Int[]
    for idx in eachindex(floes)
        isassigned(pend, idx) || continue
        id = floes.id[idx]
        for (k, mode) in pend[idx]
            nid = floes.bonds[idx][k]                   # still the neighbour id
            j   = parent_idx[nid]
            kn  = findfirst(==(id), floes.bonds[j])     # partner edge, found before clearing
            _clear_bond!(floes, idx, k, mode)
            kn === nothing || _clear_bond!(floes, j, kn, mode)
            push!(changed, idx)
            push!(changed, j)
        end
    end

    # --- Step 4. Copy updated bonds to the ghosts of the changed floes ---
    for p in unique!(changed)
        gl = id_to_indices[floes.id[p]]
        length(gl) > 1 || continue
        for g in gl
            g == p && continue
            copyto!(floes.bonds[g], floes.bonds[p])
        end
    end

    return
end