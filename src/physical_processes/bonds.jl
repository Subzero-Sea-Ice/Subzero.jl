"""
Functions needed for calculating bond forces and torques
"""

# Remove bonds to floes that no longer exist (e.g. floes that left through open boundaries)
function update_bonds!(floes)
    existing_ids = Set(floes.id)
    for (index, bonds) in enumerate(floes.bonds)
        for j in eachindex(bonds)
            if bonds[j] > 0 && !(bonds[j] in existing_ids)
                floes.bonds[index][j] = 0
            end
        end
    end
end

# Midpoint of edge k of a ring (ring[k] -> ring[k+1]), no allocation
@inline edge_mid(ring, k) = ((ring[k][1] + ring[k+1][1]) / 2, (ring[k][2] + ring[k+1][2]) / 2)

# Differences in velocity (neighbour - floe) at point `mid`, including rigid rotation
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
    domain,        # unused, kept so existing callers do not change
    consts,
    Δt,
    curr_time,     # unused, kept so existing callers do not change
    bond_settings,
) where {FT<:AbstractFloat}

    ζ_b = FT(bond_settings.ζ_b)
    ϵ   = eps(FT)
    E   = FT(consts.E)
    ν   = FT(consts.ν)
    dt  = FT(Δt)

    # --- Step 1. Remove bonds to floes that no longer exist ---
    update_bonds!(floes)

    # --- Step 2. Lookup tables, built once per timestep ---
    # id -> all indices with that id (parent + ghosts); id -> index of the parent
    id_to_indices = Dict{Int, Vector{Int}}()
    parent_idx    = Dict{Int, Int}()
    for i in eachindex(floes)
        push!(get!(id_to_indices, floes.id[i], Int[]), i)
        floes.ghost_id[i] == 0 && (parent_idx[floes.id[i]] = i)
    end

    # --- Step 3. Per-bond forces (each bond is computed once, from the floe with the smaller id) ---
    Threads.@threads for idx in eachindex(floes)
        # reset per-edge xy forces and stiffnesses (nt history is kept)
        fxy_all = floes.ind_bond_force_xy[idx]
        kst_all = floes.ind_bond_stiff_nt[idx]
        for k in eachindex(fxy_all)
            fill!(fxy_all[k], zero(FT))
            fill!(kst_all[k], zero(FT))
        end

        floes.ghost_id[idx] == 0 || continue
        id      = floes.id[idx]
        bonds_i = floes.bonds[idx]

        for k in eachindex(bonds_i)
            nid = bonds_i[k]
            (nid <= 0 || nid < id) && continue          # process each bond once: id < neighbor id

            jp = get(parent_idx, nid, 0)
            jp == 0 && continue
            kn = findfirst(==(id), floes.bonds[jp])     # edge of the neighbour that points back to us
            kn === nothing && continue

            # --- choose closest (parent/ghost of i) x (parent/ghost of j) pair ---
            best_a, best_b = idx, jp
            ia = id_to_indices[id]
            ib = id_to_indices[nid]
            if length(ia) > 1 || length(ib) > 1         # ghosts exist only near periodic boundaries
                best_d2 = typemax(FT)
                for a in ia
                    ma = edge_mid(floes.coords[a][1], k)
                    for b in ib
                        mb = edge_mid(floes.coords[b][1], kn)
                        d2 = (mb[1] - ma[1])^2 + (mb[2] - ma[2])^2
                        if d2 < best_d2
                            best_d2 = d2
                            best_a, best_b = a, b
                        end
                    end
                end
            end

            # --- bond geometry from the SELECTED floes (parent or ghost) ---
            ring_a = floes.coords[best_a][1]
            ring_b = floes.coords[best_b][1]
            a1, a2 = ring_a[k],  ring_a[k + 1]
            b1, b2 = ring_b[kn], ring_b[kn + 1]
            p1x = (a1[1] + b2[1]) / 2
            p1y = (a1[2] + b2[2]) / 2
            p2x = (a2[1] + b1[1]) / 2
            p2y = (a2[2] + b1[2]) / 2
            mx  = (p1x + p2x) / 2
            my  = (p1y + p2y) / 2
            Lb  = hypot(p2x - p1x, p2y - p1y)

            # tangent and normal (normal points from floe i towards floe j)
            inv_L = 1 / (Lb + ϵ)
            tx = (p2x - p1x) * inv_L
            ty = (p2y - p1y) * inv_L
            nx, ny = -ty, tx
            ci = floes.centroid[best_a]
            cj = floes.centroid[best_b]
            dxc = cj[1] - ci[1]
            dyc = cj[2] - ci[2]
            if nx * dxc + ny * dyc < 0
                nx, ny = -nx, -ny
            end

            # relative velocity at the bond midpoint
            v_rel_x, v_rel_y = get_vel_diff(floes, best_a, best_b, (mx, my))
            v_rel_n = v_rel_x * nx + v_rel_y * ny
            v_rel_t = v_rel_x * tx + v_rel_y * ty

            # bond properties
            h_eff  = min(floes.height[best_a], floes.height[best_b])
            area_b = Lb * h_eff
            len_b  = sqrt(dxc^2 + dyc^2) + ϵ
            I_z    = (h_eff * Lb^3) / 12
            k_b    = E * area_b / (len_b + ϵ)
            G_b    = (E / (2 * (1 + ν))) * area_b / (len_b + ϵ)
            k_θ    = E * I_z / (len_b + ϵ)

            # history (stored on the parent floe i, edge k)
            fnt     = floes.ind_bond_force_nt[idx][k]
            nrm_old = floes.ind_bond_normal_old[idx][k]
            Fn_old, Ft_old = fnt[1], fnt[2]
            nx_old, ny_old = nrm_old[1], nrm_old[2]
            if nx_old == 0 && ny_old == 0
                nx_old, ny_old = nx, ny
            end
            cos_dθ = nx_old * nx + ny_old * ny
            sin_dθ = nx_old * ny - ny_old * nx
            Fn_old_rot =  Fn_old * cos_dθ + Ft_old * sin_dθ
            Ft_old_rot = -Fn_old * sin_dθ + Ft_old * cos_dθ

            Fn_bond = k_b * v_rel_n * dt + Fn_old_rot * ζ_b
            Ft_bond = G_b * v_rel_t * dt + Ft_old_rot * ζ_b
            Fx_bond = Fn_bond * nx + Ft_bond * tx
            Fy_bond = Fn_bond * ny + Ft_bond * ty

            # store state on floe i
            nrm_old[1] = nx
            nrm_old[2] = ny
            fnt[1] = Fn_bond
            fnt[2] = Ft_bond
            fxy = fxy_all[k]
            fxy[1] = Fx_bond
            fxy[2] = Fy_bond
            kst = kst_all[k]
            kst[1] = k_b
            kst[2] = G_b
            prop = floes.ind_bond_prop[idx][k]
            prop[1] = area_b
            prop[2] = I_z
            prop[3] = Lb

            # bending moment
            Δω_z = floes.ξ[best_b] - floes.ξ[best_a]
            floes.ind_bond_bm[idx][k] = k_θ * Δω_z * dt + floes.ind_bond_bm[idx][k] * ζ_b
        end
    end

    # --- Step 4. Mirror forces onto the floe with the larger id (action = -reaction) ---
    # Reads only entries written in step 3 (smaller id), writes only its own entries: thread safe.
    Threads.@threads for idx in eachindex(floes)
        floes.ghost_id[idx] == 0 || continue
        id      = floes.id[idx]
        bonds_i = floes.bonds[idx]
        for k in eachindex(bonds_i)
            nb = bonds_i[k]
            (nb <= 0 || id < nb) && continue            # only id > neighbor id
            j  = get(parent_idx, nb, 0)
            j == 0 && continue
            kn = findfirst(==(id), floes.bonds[j])
            kn === nothing && continue

            xy_i,  xy_j  = floes.ind_bond_force_xy[idx][k], floes.ind_bond_force_xy[j][kn]
            nt_i,  nt_j  = floes.ind_bond_force_nt[idx][k], floes.ind_bond_force_nt[j][kn]
            st_i,  st_j  = floes.ind_bond_stiff_nt[idx][k], floes.ind_bond_stiff_nt[j][kn]
            pr_i,  pr_j  = floes.ind_bond_prop[idx][k],     floes.ind_bond_prop[j][kn]
            xy_i[1] = -xy_j[1];  xy_i[2] = -xy_j[2]
            nt_i[1] = -nt_j[1];  nt_i[2] = -nt_j[2]
            st_i[1] =  st_j[1];  st_i[2] =  st_j[2]
            pr_i[1] =  pr_j[1];  pr_i[2] =  pr_j[2];  pr_i[3] = pr_j[3]
            floes.ind_bond_bm[idx][k] = -floes.ind_bond_bm[j][kn]
        end
    end

    # --- Step 5. Total bond force, torque and stress for every floe (one pass, no allocations) ---
    Threads.@threads for idx in eachindex(floes)
        sx = zero(FT); sy = zero(FT); trq = zero(FT)
        sxx = zero(FT); sxy = zero(FT); syy = zero(FT)
        parent = floes.ghost_id[idx] == 0
        ring   = floes.coords[idx][1]
        cx, cy = floes.centroid[idx][1], floes.centroid[idx][2]
        fxy_all = floes.ind_bond_force_xy[idx]
        bm_all  = floes.ind_bond_bm[idx]
        bonds_all = floes.bonds[idx]

        for k in eachindex(fxy_all)
            Fx = fxy_all[k][1]
            Fy = fxy_all[k][2]
            sx += Fx
            sy += Fy
            if parent
                mx, my = edge_mid(ring, k)
                rx, ry = mx - cx, my - cy
                trq += (rx * Fy - ry * Fx) + (bonds_all[k] > 0 ? bm_all[k] : zero(FT))
                sxx += Fx * rx
                sxy += Fx * ry + Fy * rx
                syy += Fy * ry
            end
        end

        floes.bond_force[idx][1] = sx
        floes.bond_force[idx][2] = sy
        floes.bond_trq[idx] = trq

        inv_v = 1 / (floes.area[idx] * floes.height[idx] + ϵ)
        S = floes.stress_bond[idx]
        S[1, 1] = sxx * inv_v
        S[1, 2] = FT(0.5) * sxy * inv_v
        S[2, 1] = S[1, 2]
        S[2, 2] = syy * inv_v
    end

    return
end