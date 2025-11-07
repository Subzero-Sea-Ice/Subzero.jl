"""
Functions needed for calculating bond forces and torques
"""
function update_bonds!(floes)
    existing_ids = Set(floes.id)
    #println("existing_ids::", existing_ids )
    for (index, bonds) in enumerate(floes.hex_bonds)
        for j in eachindex(bonds)
            if bonds[j] != 0 && !(bonds[j] in existing_ids)
                floes.hex_bonds[index][j] = 0
            end
        end
    end
end

function polygon_edges(poly_coords::Vector{Vector{Float64}})
    n = length(poly_coords)
    pts = [GB.Point2(c[1], c[2]) for c in poly_coords]
    [GB.Line(pts[i], pts[mod1(i+1, n)]) for i in 1:n]
end

function midpoint(edge::GB.Line)
    p1, p2 = edge.points
    GB.Point2((p1[1] + p2[1])/2, (p1[2] + p2[2])/2)
end

# Function to check if point is inside polygon (outer ring)
function point_in_polygon(p::GB.Point2, poly_coords::Vector{Vector{Float64}})
    x, y = p[1], p[2]
    n = length(poly_coords)
    inside = false
    j = n
    tol = 1e-9

    for i in 1:n
        xi, yi = poly_coords[i]
        xj, yj = poly_coords[j]

        # --- Check if point lies on this edge ---
        cross = (x - xi) * (yj - yi) - (y - yi) * (xj - xi)
        if abs(cross) < tol
            if (min(xi, xj) - tol ≤ x ≤ max(xi, xj) + tol) &&
               (min(yi, yj) - tol ≤ y ≤ max(yi, yj) + tol)
                return true  # point is on the boundary
            end
        end

        # --- Ray casting for inside/outside test ---
        intersect = ((yi > y) != (yj > y)) &&
                    (x < (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)

        inside ⊻= intersect
        j = i
    end

    return inside
end

#Function to calculate differences in velocities at the endpoints of a bond
function get_vel_diff(floes, idx_floe, idx_neighbour, c_floe, c_neighbour)
    u_floe = floes.u[idx_floe] + floes.ξ[idx_floe] * (c_floe[1] - floes.centroid[idx_floe][1])
    v_floe = floes.v[idx_floe] + floes.ξ[idx_floe] * (c_floe[2] - floes.centroid[idx_floe][2])
    u_neighbour = floes.u[idx_neighbour] + floes.ξ[idx_neighbour] * (c_neighbour[1] - floes.centroid[idx_neighbour][1])
    v_neighbour = floes.v[idx_neighbour] + floes.ξ[idx_neighbour] * (c_neighbour[2] - floes.centroid[idx_neighbour][2])
    Δu = u_neighbour - u_floe
    Δv = v_neighbour - v_floe
    return Δu, Δv
end

function timestep_hexagonal_bonds!(
    floes::StructArray{<:Floe{FT}},
    domain,
    consts,
    curr_time,
    hexagonal_bond_settings,
) where {FT<:AbstractFloat}
    nfaces = hexagonal_bond_settings.nfaces
    k_ice  = hexagonal_bond_settings.elastic_constant
    c_ice  = hexagonal_bond_settings.damp_coeff
    update_bonds!(floes)
    for id in eachindex(floes)
        for i in eachindex(floes.ind_bond_force[id])
            floes.ind_bond_force[id][i] .= FT(0)  # set both Fx and Fy to zero
        end
    end
    for idx in eachindex(floes)
        id = floes.id[idx]
        if floes.ghost_id[idx] == 0
            fill!(floes.bond_force[idx], FT(0))  # Reset both Fx and Fy to zero
            i_indices = findall(==(id), floes.id)  # actual + ghost indices
            edges_i = polygon_edges(floes.coords[i_indices[1]][1])
            num_edges_i = length(edges_i) - 1  
            for edge_idx in 1:num_edges_i
                neighbor = floes.hex_bonds[i_indices[1]][edge_idx]  
                Fx_bond, Fy_bond = 0.0, 0.0
                if neighbor != 0 && id < neighbor
                    neighbor_indices = findall(==(neighbor), floes.id)
                    best_dx, best_dy = 0.0, 0.0
                    best_dist = Inf
                    best_pair = nothing
                    for i_idx in i_indices
                        edges_i2 = polygon_edges(floes.coords[i_idx][1])
                        edge_i = edges_i2[edge_idx]
                        mid_i = midpoint(edge_i)
                        for neighbor_idx in neighbor_indices
                            edges_neighbor = polygon_edges(floes.coords[neighbor_idx][1])
                            neighbor_edge_idx = findfirst(x -> x == id, floes.hex_bonds[neighbor_idx])
                            edge_neighbor = edges_neighbor[neighbor_edge_idx]
                            mid_neighbor = midpoint(edge_neighbor)
                            dx = mid_neighbor[1] - mid_i[1]
                            dy = mid_neighbor[2] - mid_i[2]
                            dist = sqrt(dx^2 + dy^2)
                            in_neighbor = point_in_polygon(mid_i, floes.coords[neighbor_idx][1])
                            in_i = point_in_polygon(mid_neighbor, floes.coords[i_idx][1])
                            if in_neighbor || in_i
                                dx, dy = 0.0, 0.0
                                dist = 0.0
                            end
                            if dist < best_dist
                                best_dx, best_dy = dx, dy
                                best_dist = dist
                                best_pair = (i_idx, edge_idx, neighbor_idx, neighbor_edge_idx, mid_i, mid_neighbor)
                            end
                        end
                    end
                    # assign forces once for the best pair
                    if best_pair !== nothing
                        i_idx2, edge_idx2, neighbor_idx, neighbor_edge_idx, mid_i2, mid_neighbor2 = best_pair
                        floe_idx_act = findfirst(x -> x == floes.id[i_idx2], floes.id)
                        Δu, Δv = get_vel_diff(floes, i_idx2, neighbor_idx, mid_i2, mid_neighbor2)
                        Fx_bond = best_dx * k_ice +  c_ice * Δu
                        Fy_bond = best_dy * k_ice +  c_ice * Δv
                        floes.ind_bond_force[floe_idx_act][edge_idx2][1] = Fx_bond
                        floes.ind_bond_force[floe_idx_act][edge_idx2][2] = Fy_bond
                    end
                end 
                if neighbor != 0
                    floe_idx_act = findfirst(x -> x == id, floes.id)
                    nei_idx_act = findfirst(x -> x == neighbor, floes.id)
                    nei_edge_idx = findfirst(x -> x == id, floes.hex_bonds[nei_idx_act])
                    floes.ind_bond_force[nei_idx_act][nei_edge_idx][1] = - floes.ind_bond_force[floe_idx_act][edge_idx][1] 
                    floes.ind_bond_force[nei_idx_act][nei_edge_idx][2] = - floes.ind_bond_force[floe_idx_act][edge_idx][2] 
                end
            end    
            floes.bond_force[idx][1] = sum(first.(floes.ind_bond_force[idx]))
            floes.bond_force[idx][2] = sum(last.(floes.ind_bond_force[idx]))
            #if(curr_time == 5000)
            #    step_t = curr_time/20
            #    println(">>Time:", step_t, " >>Floe:", id,  " >>Forces:" ,floes.ind_bond_force[idx])
            #    println(">>Time:", step_t, " >>Floe:",  id,  " >>Forces:" ,floes.bond_force[idx])
            #end
        end
    end


    # Loop for calculating bond torques, stresses
    for idx in eachindex(floes)
        floes.bond_trq[idx] = FT(0.0)
        fill!(floes.stress_bond[idx], FT(0))  # Reset all stresses to zero
        floe_stress_bond =  zeros(FT, 2, 2)
        if floes.ghost_id[idx] == 0
            edges_i = polygon_edges(floes.coords[idx][1])
            num_edges_i = length(edges_i) - 1  
            for edge_idx in 1:num_edges_i
                torque_xy = 0.0
                edge_i = edges_i[edge_idx]
                mid_i = midpoint(edge_i)
                r_x = mid_i[1] - floes.centroid[idx][1]
                r_y = mid_i[2] - floes.centroid[idx][2]
                Fx_bond, Fy_bond = floes.ind_bond_force[idx][edge_idx]
                torque_xy = r_x * Fy_bond - r_y * Fx_bond
                # Accumulaye all bond torques for each floe
                floes.bond_trq[idx] += torque_xy
                # Accumulate stresses
                floe_stress_bond[1,1] += Fx_bond * r_x
                floe_stress_bond[1,2] += Fx_bond * r_y + Fy_bond * r_x
                floe_stress_bond[2,2] += Fy_bond * r_y
                #bond_id = floes.hex_bonds[idx][edge_idx]  
                #if(curr_time == 5000)
                #    step_t = curr_time/20
                #    println("Time: $step_t, Floe: $idx, Neighbor: $bond_id, Forces: $Fx_bond, $Fy_bond, Torque: $torque_xy, Total:", floes.bond_trq[idx])
                #    println("Time: $step_t, Floe: $idx, coords: ", floes.coords[idx][1])
                #    println("Time: $step_t, Floe: $idx, Neighbor: $bond_id, x,y: $mid_i")
                #    println("Time: $step_t, Floe: $idx, Neighbor: $bond_id", "<cx:>", floes.centroid[idx][1], "<cy:>", floes.centroid[idx][2])
                #    println("Time: $step_t, Floe: $idx, Neighbor: $bond_id, r_x: $r_x, r_y: $r_y")
                #end 
            end
        end
        floe_stress_bond[1, 2] *= FT(0.5)
        floe_stress_bond[2, 1] = floe_stress_bond[1, 2]
        floe_stress_bond .*= 1/(floes.area[idx] * floes.height[idx])
        floes.stress_bond[idx] .= floe_stress_bond
    end

    return
end   
