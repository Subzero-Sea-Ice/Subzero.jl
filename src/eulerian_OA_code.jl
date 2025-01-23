"""
Code for finding OA forcings using Eulerian Grid calculations. 
Might be used later for allowing users the choice of using monte carlo
or Eulerian. Not currently used. 
"""

"""
floe_grid_bounds(g, p, rmax)

Finds the bounding indices of a grid line such that a point plus and minus a maximum radius are within those grid lines.
Inputs:
    g <Vector{Real}> grid lines
    p <Real> point value
    rmax <Real> radial buffer around p
Outputs:
    min_idx: index of g such that the value of (p - rmax) is between 
             indices g[min_idx] and g[min_idx + 1]. If p is less than all of the values in g, this will be 0.
    max_idx: index of g such that the value of (p + rmax) is between 
             indices g[max_idx - 1] and g[man_idx]. If p is greater than all of the values in g, this will be the length of g.
Note: If radius is negative this will switch the minimum and maximum indices.
"""
function floe_grid_bounds(g, p, rmax)
    Δ = g[2] - g[1]
    pmin = p - rmax > g[1] ? p - rmax : g[1]
    pmax = p + rmax < g[end] ? p + rmax : g[end]
    min_idx = findmin(abs.(g .- pmin))[2]
    min_val = g[min_idx]
    if pmin < min_val
        min_idx -= 1
        min_val -= Δ
    end
    max_idx = Int(cld(pmax-pmin, Δ)) + min_idx
    if pmax > g[max_idx]
        max_idx += 1
    end
    return min_idx, max_idx
end

"""
cell_area_ratio(cell_poly, floe_poly)

Calculates the percentage of a grid square filled with a given floe.
Inputs:
    cell_poly <LibGEOS.Polygon>
    floe_poly <LibGEOS.Polygon>
Outputs:
    cell area ratio <Float> ratio of cell area filled with given floe
"""
function cell_area_ratio(cell_poly, floe_poly)
    floe_in_cell = LG.intersection(floe_poly, cell_poly)
    return GO.area(floe_in_cell)/GO.area(cell_poly)
end

"""
domain_coords(domain::Domain)
Inputs:
    domain<Domain>
Output:
    RingVec coordinates for edges of rectangular domain based off of boundary values
"""
function cell_coords(xmin, xmax, ymin, ymax)
    return [[[xmin, ymax], [xmin, ymin],
            [xmax, ymin], [xmax, ymax],
            [xmin, ymax]]]
end

"""
floe_area_ratio(floe, xg, yg)

Calculates the cell area ratio of grid squares surrounding given floe and the indicies of those grid squares within the grid defined by gridlines xg and yg.
Inputs:
    floe    <Floe>
    xg      <Vector{Float}> x grid lines
    yg      <Vector{Float}> y grid lines
Outputs:
    area_ratios <Vector{Float}> vector of area ratios for each grid cell 
                                within floe grid bounds
    xidx <Vector{Int}> x indices of grid cells - order matches area_ratios
    yidx <Vector{Int}> y indices of grid cells - order matches area_ratios
    idx <Vector{(Int, Int)}> cartesian point defining one grid cell in grid
                             - order matches area_ratios
"""
function floe_area_ratio(floe, xg, yg, t::Type{T} = Float64) where T
    floe_poly = LG.Polygon(floe.coords)
    xmin_idx, xmax_idx = floe_grid_bounds(xg, floe.centroid[1], floe.rmax)
    ymin_idx, ymax_idx = floe_grid_bounds(yg, floe.centroid[2], floe.rmax)
    nx = xmax_idx - xmin_idx
    ny = ymax_idx - ymin_idx
    area_ratios = T[]
    xidx = Int[]
    yidx = Int[]
    for i = xmin_idx:(xmax_idx-1)
        for j = ymin_idx:(ymax_idx-1)
            cell_poly = LG.Polygon(cell_coords(xg[i], xg[i+1], yg[j], yg[j+1]))
            ratio = cell_area_ratio(cell_poly, floe_poly)
            if ratio > 0.0
                push!(area_ratios, ratio)
                push!(xidx, i)
                push!(yidx, j)
            end
        end
    end
    # y values are rows and x values are columns
    idx = CartesianIndex.(Tuple.(eachrow(hcat(yidx,xidx))))
    return area_ratios, xidx, yidx, idx
end

"""
calc_OA_forcings!(m, i)

Calculate the effects on the ocean and atmpshere on floe i within the given model
and the effects of the ice floe on the ocean.

Inputs:
    floe    <Floe> floe
    m       <Model> given model
Outputs:
    None. Both floe and ocean fields are updated in-place.
Note: For floes that are completly out of the Grid, simulation will error. 
"""
function floe_OA_forcings!(floe, m, c)
    Δx = m.grid.xg[2] - m.grid.xg[1]
    Δy = m.grid.yg[2] - m.grid.yg[1]

    # Grid squares under ice floe and ice area per cell
    ma_ratio = floe.mass/floe.area
    area_ratios, xidx, yidx, idx = floe_area_ratio(floe, m.grid.xg, m.grid.yg)
    areas = area_ratios * (Δx * Δy)

    # Floe heatflux
    floe.hflx = mean(m.ocean.hflx[idx])

    # Ice velocity within each grid square
    lx = m.grid.xc[xidx] .- floe.centroid[1]
    ly = m.grid.yc[yidx] .- floe.centroid[2]
    uice = floe.u .- ly*floe.ξ
    vice = floe.v .+ lx*floe.ξ

    # Force on ice from atmopshere
    uatm = m.wind.u[idx]
    vatm = m.wind.v[idx]
    fx_atm = (c.ρa * c.Cd_ia * sqrt.(uatm.^2 + vatm.^2) .* uatm) .* areas
    fy_atm = (c.ρa * c.Cd_ia * sqrt.(uatm.^2 + vatm.^2) .* vatm) .* areas

    # Force on ice from pressure gradient
    fx_pressure∇ = -ma_ratio * c.f .* m.ocean.v[idx] .* areas
    fy_pressure∇ = ma_ratio * c.f .* m.ocean.u[idx] .* areas

    # Force on ice from ocean
    Δu_OI = m.ocean.u[idx] .- uice
    Δv_OI = m.ocean.v[idx] .- vice
    τx_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (cos(c.turnθ) .* Δu_OI .- sin(c.turnθ) * Δv_OI)
    τy_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (sin(c.turnθ) .* Δu_OI .+ cos(c.turnθ) * Δv_OI)
    fx_ocn = τx_ocn .* areas
    fy_ocn = τy_ocn .* areas

    # Sum above forces and find torque
    fx = fx_atm .+ fx_pressure∇ .+ fx_ocn
    fy = fy_atm .+ fy_pressure∇ .+ fy_ocn
    trq = lx.*fy .- ly.*fx  # are these signs right?

    # Add coriolis force to total foces
    fx .+= ma_ratio * c.f * floe.v * areas
    fy .-= ma_ratio * c.f * floe.u * areas

    # Sum forces on ice floe
    floe.fxOA = sum(fx)
    floe.fyOA = sum(fy)
    floe.torqueOA = sum(trq)

    # TODO: Not thread safe
    # Update ocean stress fields with ice on ocean stress
    m.ocean.τx[idx] .= m.ocean.τx[idx].*(1 .- area_ratios) .- τx_ocn.*area_ratios
    m.ocean.τy[idx] .= m.ocean.τy[idx].*(1 .- area_ratios) .- τy_ocn.*area_ratios

    # Update sea-ice fraction
    m.ocean.si_frac[idx] .+= area_ratios
    return
end

# @testset "Simulation Creation" begin
#     @testset "Floe Area Ratio" begin
#         @test Subzero.floe_grid_bounds(collect(0:2:10), 5, 1.25) == (2, 5)
#         @test Subzero.floe_grid_bounds(collect(0:2:10), 5, 2.25) == (2, 5)
#         @test Subzero.floe_grid_bounds(collect(0:2:10), 0.75, 1.25) == (1, 2)
#         @test Subzero.floe_grid_bounds(collect(0:2:10), 9.25, 1.25) == (5, 6)

#         cell_poly = LibGEOS.Polygon([[[0., 10.], [0., 0.], [10., 0.],
#                                       [10., 10.], [0., 10.]]])
#         floe_poly1 = LibGEOS.Polygon([[[2., 8.], [2., 2.], [8., 2.],
#                                          [8., 8.], [2., 8.]]])
#         floe_poly2 = LibGEOS.Polygon([[[9., 5.], [9., 2.], [12., 2.],
#                                          [12., 5.], [9., 5.]]])
#         @test Subzero.cell_area_ratio(cell_poly, floe_poly1) == 0.36
#         @test Subzero.cell_area_ratio(cell_poly, floe_poly2) == 0.03
#         @test Subzero.cell_area_ratio(cell_poly,
#                 Subzero.translate(floe_poly2, [1.0, 0.0])) == 0.0
#         floe_poly3 = LibGEOS.Polygon([[[3., 8.], [3., 3.], [8., 3.],
#                                        [8., 8.],[3., 8.]]])
#         floe3 = Subzero.Floe(floe_poly3, 0.25, 0.0)
#         area_ratio3, _, _, idx1 = Subzero.floe_area_ratio(floe3, collect(-5.:5.:10.), collect(0.:5.:10.))
#         @test area_ratio3[1] == 4/25
#         @test area_ratio3[2] == 6/25
#         @test area_ratio3[3] == 6/25
#         @test area_ratio3[4] == 9/25

#         floe_poly4 = LibGEOS.Polygon([[[2., 2.], [8., 2.], [8., 8.], [2., 2.]]])
#         floe4 = Subzero.Floe(floe_poly4, 0.25, 0.0)
#         area_ratio4, _, _, idx1 = Subzero.floe_area_ratio(floe4, collect(-5.:5.:10.), collect(0.:5.:10.))
#         @test area_ratio4[1] == 0.18
#         @test area_ratio4[2] == 0.36
#         @test area_ratio4[3] == 0.18
#     end
# end