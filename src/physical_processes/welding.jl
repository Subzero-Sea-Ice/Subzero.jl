"""
    bin_floe_centroids(floes, grid, domain, Nx, Ny)

Split floe locations into a grid of Nx by Ny by floe centroid location
Inputs:
    floes       <StructArray{Floe}> simulation's list of floes
    grid        <RegRectilinearGrid> simulation's grid
    domain      <Domain> simulation's domain
    Nx          <Int> number of grid cells in the x-direction to split domain
                    into for welding groups
    Ny          <Int> number of grid cells in the y-direction to split domain
                    into for welding groups
Outputs:
    floe_bins   <Matrix{Vector{Int}}> Nx by Ny matrix where each element is a
                    list of floe indices whose centroids are in the
                    corresponding section of the grid. May also have elements of
                    value 0 if there are less than average number of floes in 
                    the section.
    nfloes      <Maxtrix{Int}> Nx by Ny matrix where each element is the total
                    number of indices within floe_bins[Nx, Ny] that are
                    non-zeros and represent a floe within the grid section.
"""
function bin_floe_centroids(floes, grid, domain, Nx, Ny)
    @assert Nx > 0 && Ny > 0 "Can't bin centroids without bins."
    # Find average number of floes per bin if floes were spread evenly
    floes_per_bin = ceil(Int, length(floes) / (Nx * Ny))
    # Create bins with vectors of zeros
    floe_bins = [zeros(Int, floes_per_bin) for _ in 1:Nx, _ in 1:Ny]
    nfloes = zeros(Int, Nx, Ny)
    # Determine length of new grid bins
    Lx = grid.xf - grid.x0
    Ly = grid.yf - grid.y0
    Δx = Lx / Nx
    Δy = Ly / Ny
    for i in eachindex(floes)
        # Determine if centroid is within bounds
        xp, yp = floes.centroid[i]
        !in_bounds(xp, yp, grid, domain.north, domain.east) && break
        # Determine which grid bin centroid is in and shift based on boundary
        xidx = grid_cell_index(xp, Δx, grid.x0)
        xidx = (xp ≤ grid.x0) ? 1 : xidx
        xidx = (xp ≥ grid.xf) ? Nx : xidx
        yidx = grid_cell_index(yp, Δy, grid.y0)
        yidx = (yp ≤ grid.y0) ? 1 : yidx
        yidx = (yp ≥ grid.yf) ? Ny : yidx
        # Add floe to bin list
        nfloes[xidx, yidx] += 1
        if nfloes[xidx, yidx] > floes_per_bin
            push!(floe_bins[xidx, yidx], i)
        else
            floe_bins[xidx, yidx][nfloes[xidx, yidx]] = i
        end
    end
    return floe_bins, nfloes
end

"""
    timestep_welding!(
        floes,
        max_floe_id,
        grid,
        domain,
        Nx,
        Ny,
        weld_settings::WeldSettings{FT},
        floe_settings
        Δt,
        rng,
    )

Weld floes within sections of the domain that meet overlap and size criteria
together, ensuring resulting floe doesn't surpass maximum floe area. 
Inputs:
    floes               <StructArray{Floe}> simulation's list of floes
    max_floe_id         <Int> maximum floe ID before this welding
    grid                <RegRectilinearGrid> simulation's grid
    domain              <Domain> simulation's domain
    Nx                  <Int> number of grid cells in the x-direction to split
                            domain into for welding groups
    Ny                  <Int> number of grid cells in the y-direction to split
                            domain into for welding groups
    weld_settings       <WeldSettings> welding settings
    floe_settings       <FloeSettings> sim's settings for making new floes
    consts              <Consts> simulation's constants
    Δt                  <Int> length of timestep in seconds
    rng                 <RandomNumberGenerator> simulation's rng
Outputs:
    Returns nothing. Welds groups of floes together that meet requirments. Floes
    that are fused into other floes are marked for removal.
"""
function timestep_welding!(
    floes,
    max_floe_id,
    grid,
    domain,
    weld_settings::WeldSettings{FT},
    floe_settings,
    weld_idx,
    Δt,
    rng = Xoshiro(),
) where FT <: AbstractFloat
    # Seperate floes into groups based on centroid location in grid
    Nx, Ny = weld_settings.Nxs[weld_idx], weld_settings.Nys[weld_idx]
    floe_bins, floes_per_bin = bin_floe_centroids(floes, grid, domain, Nx, Ny)
    for k in eachindex(floe_bins)  # should be able to multithread
        bin = floe_bins[k]
        nfloes = floes_per_bin[k]
        weld_group = [(0, FT(0)) for _ in 1:nfloes]
        @views for i in bin[1:nfloes]
            # Reset the welding group
            weld_slot_idx = 0
            #=
            Determine which floes can weld with first floe in bin given it isn't
            too large and is still active in the simulation
            =#
            if (floes.status[i].tag == active &&
                floes.area[i] < weld_settings.max_weld_area
            )
                for j in bin[1:nfloes]
                    if (
                        # don't re-check pairs, but check new floes from fusion
                        i != j &&
                        (i < j || floes.id[j] > max_floe_id) &&
                        floes.status[i].tag == active &&
                        floes.status[j].tag == active &&
                        floes.area[i] < weld_settings.max_weld_area &&
                        floes.area[j] < weld_settings.max_weld_area &&
                        potential_interaction(  # floes must be interacting
                            floes.centroid[i], floes.centroid[j],
                            floes.rmax[i], floes.rmax[j]
                        )
                    )
                        # Find intersection area
                        inter_area = FT(sum(GO.area, intersect_polys(floes.poly[i], floes.poly[j], FT); init = 0.0))
                        # Probability two floes will weld
                        weld_prob = weld_settings.welding_coeff *
                            (inter_area / floes.area[i])
                        # Area of floes welded together
                        union_area = floes.area[i] + floes.area[j] - inter_area
                        if (
                            inter_area > 0 &&  # must overlap
                            weld_prob > rand(rng) &&  # meets weld probability
                            weld_settings.min_weld_area < union_area &&
                            weld_settings.max_weld_area > union_area
                        )
                            # Add floe j to list of floes i can weld with
                            weld_slot_idx += 1
                            weld_group[weld_slot_idx] = (j, inter_area)
                        end
                    end
                end
            end
            # Sort floes that floe i will weld with by intersection area
            sort!(weld_group[1:weld_slot_idx], by = x  -> last(x), rev = true)
            for idx in 1:weld_slot_idx
                j, ij_inter_area = weld_group[idx]
                # If new floe i will be too large, stop welding
                new_area = floes.area[i] + floes.area[j] - ij_inter_area
                new_area > weld_settings.max_weld_area && break
                # Weld floe i and j, replacing floe i with welded floe
                fuse_two_floes!(
                    get_floe(floes, i),
                    get_floe(floes, j),
                    Δt,
                    floe_settings,
                    max_floe_id,
                    rng,
                )
                # can't update to real value due to multi-threading - do below
                floes.id[i] = -1
            end
        end
    end
    # Update ids of floes outside of multithreading loop
    for i in eachindex(floes)
        if floes.id[i] == -1
            max_floe_id += 1
            floes.id[i] = max_floe_id
        end
    end
    return max_floe_id
end