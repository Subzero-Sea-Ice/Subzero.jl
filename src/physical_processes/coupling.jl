# Functions needed for coupling between floes, ocean, and atmosphere

#=
Find the index of given point's cartesian value (in either the x or y direction)
within the simulation grid. 
=#
grid_cell_index(p, Δg, g0) = floor(Int, (p - g0)/Δg) + 1

#=
Find index of the grid cell of the given RegRectilinearGrid that the given
x-coordinate and y-coordinate falls within.
Method depends on grid being a regular rectilinear grid.
Inputs:
    xp      <AbstractFloat> x-coordinates of point
    yp      <AbstractFloat> y-coordinate of point
    grid    <RegRectilinearGrid> simulation grid
Outputs:
    xidx    <AbstractFloat> x-index of grid cell x-point is within - this is the
                column
    yidx    <AbstractFloat> y-index of grid cell y-point is within - this is the
                row
Note:
    Points can be outside of the grid, so index can be less than 1 or greater
    than the number of grid cells
=#
function find_grid_cell_index(xp, yp, grid::RegRectilinearGrid)
    xidx = floor(Int, (xp - grid.x0) / grid.Δx) + 1
    yidx = floor(Int, (yp - grid.y0) / grid.Δy) + 1
    return xidx, yidx
end


#=
Find index of the cell centered on grid lines of the given RegRectilinearGrid
that the given x-coordinate and y-coordinate falls within.
This cell is centered around the grid lines, so it is a shifted grid cell
by half a cell. Method depends on grid being a regular rectilinear grid.
Inputs:
    xp      <AbstractFloat> x-coordinates of point
    yp      <AbstractFloat> y-coordinate of point
    grid    <RegRectilinearGrid> simulation grid
Outputs:
    xidx    <AbstractFloat> x-index of grid cell (centered on grid lines) x-point
                is within - this is the column
    yidx    <AbstractFloat> y-index of grid cell (centered on grid lines) y-point
                is within - this is the row
Note:
    Points can be outside of the grid, so index can be less than 1 or greater
    than the number of grid lines in a given direction.
=#
function find_center_cell_index(xp, yp, grid::Union{RegRectilinearGrid, dev_RegRectilinearGrid})
    xidx = floor(Int, (xp - grid.x0)/(grid.Δx) + 0.5) + 1
    yidx = floor(Int, (yp - grid.y0)/(grid.Δy) + 0.5) + 1
    return xidx, yidx
end

#=
With all non-periodic boundaries, points outside of the grid in both the x and y
are defined to be out of bounds since these points can't be interpolated as we
don't have any information on the ocean outside of the grid.

Returns a Boolean that is true if both xr and yr are within domain boundaries 
and false otherwise.
=#
function in_bounds(
    xr,
    yr,
    grid,
    ::NonPeriodicBoundary,
    ::NonPeriodicBoundary,
)
    return (grid.x0 <= xr <= grid.xf) && 
        (grid.y0 <= yr <= grid.yf)
end

#=
With the north/south non-periodic boundaries, points outside of the grid in the
y-direction are defined to be out of bounds since these points can't be
interpolated as we don't have any information on the ocean outside of the grid.
Returns a Boolean that is true if yr is within domain boundaries, and false otherwise.
=#
function in_bounds(
    xr,
    yr,
    grid,
    ::NonPeriodicBoundary,
    ::PeriodicBoundary,
)
    return grid.y0 <= yr <= grid.yf
end

#=
With the east/west non-periodic boundaries, points outside of the grid in the
x-direction are defined to be out of bounds since these points can't be
interpolated as we don't have any information on the ocean outside of the grid.
Returns a Boolean that is true if xr is within domain boundaries, and false otherwise.
=#
function in_bounds(
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::NonPeriodicBoundary,
)
    return grid.x0 <= xr <= grid.xf
end

#=
With all periodic boundaries, all points are considered to be in-bounds.
Returns a Boolean that is true regardless of point values.
=#
function in_bounds(
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::PeriodicBoundary,
)
    return true
end

#-------------- Interpolation of Ocean and Atmosphere --------------#

"""
    axis_stencil()

Compute the indexes of the closest grid lines to a given point within a floe

## _Positional arguments_
- `x`: coordinate of the floe point
- `x0`: coordinate of the first grid line on this axis
- `Δx`: distance between neighbouring grid lines
- `N`: number of grid cells on an axis
- boundary type

## _Returns_
- `mod(i, N) + 1`: the index of the line before the point
- `mod(i + 1, N) + 1`: the index of the line after the point
- `s - k`: the weight within [0,1[
"""
function axis_stencil(x, x0, Δx, N, ::PeriodicBoundary)
    s = (x - x0) / Δx #point distance from firs line in cells 
    k = floor(s)      #number of cells passed
    i = unsafe_trunc(Int, k)   # k as an integer, with no error handling
    return mod(i, N) + 1, mod(i + 1, N) + 1, s - k
end

#=
Non-periodic boundary: the field has N + 1 grid lines and no wrap-around.
Points are already filtered by `in_bounds`, so x is within [x0, x0 + N * Δx].
The cell index is clamped to [0, N - 1] so that a point exactly on the last grid
line uses lines N and N + 1 with weight 1, instead of reading past the array.
=#
function axis_stencil(x, x0, Δx, N, ::NonPeriodicBoundary)
    s = (x - x0) / Δx
    i = clamp(unsafe_trunc(Int, floor(s)), 0, N - 1)
    return i + 1, i + 2, s - i
end

# Compute the value of one or more fields at a point (x, y) using a bilinear interpolation
# Inputs:
# - x, y = point coordinates
# - grid = model's grid on the device
# - domain = model's domain on the device
# - fields = a tuple of matrices such as (atmos.u, atmos.v)
# Returns: 
# - a tuple of numbers, one per matrix, in the same order as fields.
function interpolate_fields(x, y, grid, domain, fields::Tuple)
    i1, i2, tx = axis_stencil(x, grid.x0, grid.Δx, grid.Nx, domain.east)
    j1, j2, ty = axis_stencil(y, grid.y0, grid.Δy, grid.Ny, domain.north)
    return map(fields) do f
        (1 - ty) * ((1 - tx) * f[i1, j1] + tx * f[i2, j1]) +
             ty  * ((1 - tx) * f[i1, j2] + tx * f[i2, j2])
    end
end

# Compute the stress (force per unit area) applied by the wind to the floe at one point.
# Inputs:
# - x,y = point coordinates
# - upoint, vpoint = floe's velocity at this point in x and y, including it's rotation
# - atmos = atmosphere on the device
# - grid, domain on the device
# - constants for air density and ice-atmosphere drag coefficient
# Returns:
# - stress on the floe point on x and y
function atmosphere_forcing(x, y, upoint, vpoint, atmos, grid, domain, c)
    uatm, vatm = interpolate_fields(x, y, grid, domain, (atmos.u, atmos.v))
    Δu_AI = uatm - upoint
    Δv_AI = vatm - vpoint
    norm = sqrt(Δu_AI^2 + Δv_AI^2)
    τx_atm = c.ρa * c.Cd_ia * norm * Δu_AI
    τy_atm = c.ρa * c.Cd_ia * norm * Δv_AI
    return τx_atm, τy_atm
end

# Compute ocean's stress and read the heat flux on a floe point
# Input:
# - x, y = point coordinates
# - upoint, vpoint = floe's velocity at this point in x and y, including it's rotation
# - ocean = ocean on the device
# - ma_ratio = floe mass divided by its area
# - grid, domain on the device
# - constants for ocean density and ice-ocean drag coefficient, ocean turning angle, and coriolis frequency
# Returns:
# - drag stress from the ocean along x and y
# - the pressure gradient stress
# - heat flux at the floe point
function ocean_forcing(x, y, upoint, vpoint, ocean, ma_ratio, grid, domain, c)
    uocn, vocn, hflx_factor = interpolate_fields(x, y, grid, domain,
                                                 (ocean.u, ocean.v, ocean.hflx_factor))
    Δu_OI = uocn - upoint
    Δv_OI = vocn - vpoint
    norm = sqrt(Δu_OI^2 + Δv_OI^2)
    τx_ocn = c.ρo * c.Cd_io * norm * (cos(c.turnθ) * Δu_OI - sin(c.turnθ) * Δv_OI)
    τy_ocn = c.ρo * c.Cd_io * norm * (sin(c.turnθ) * Δu_OI + cos(c.turnθ) * Δv_OI)
    τx_pressure∇ = -ma_ratio * c.f * vocn
    τy_pressure∇ = ma_ratio * c.f * uocn
    return τx_ocn, τy_ocn, τx_pressure∇, τy_pressure∇, hflx_factor
end

#-------------- Effects of Ice and Atmosphere on Ocean --------------#
#=
    check_cell_bounds(..., ::PeriodicBoundary, ::PeriodicBoundary)

Return cell bounding values as is given the domain is doubley periodic and thus
the cell can extend beyond the grid as it will simply wrap back around into grid
through opposite periodic boundary.
=#
function check_cell_bounds(
    xmin,
    xmax,
    ymin,
    ymax,
    grid,
    ::PeriodicBoundary,
    ::PeriodicBoundary,
)
    return xmin, xmax, ymin, ymax
end

#=
    check_cell_bounds(..., ::NonPeriodicBoundary, ::PeriodicBoundary)

Trim cell bound in the north-south direction if it exends past grid due to
non-periodic boundary pair. 
=#
function check_cell_bounds(
    xmin,
    xmax,
    ymin,
    ymax,
    grid,
    ::NonPeriodicBoundary,
    ::PeriodicBoundary,
)
    ymin = ymin < grid.y0 ?
        grid.y0 :
        (ymin > grid.yf ? grid.yf : ymin)

    ymax = ymax > grid.yf ?
        grid.yf :
        (ymax < grid.y0 ? grid.y0 : ymax)
    return xmin, xmax, ymin, ymax
end

#=
    check_cell_bounds(..., ::PeriodicBoundary, ::NonPeriodicBoundary)

Trim cell bound in the east-west direction if it exends past grid due to
non-periodic boundary pair.
=#
function check_cell_bounds(
    xmin,
    xmax,
    ymin,
    ymax,
    grid,
    ::PeriodicBoundary,
    ::NonPeriodicBoundary,
)
    xmin = xmin < grid.x0 ?
        grid.x0 :
        (xmin > grid.xf ? grid.xf : xmin)

    xmax = xmax > grid.xf ?
        grid.xf :
        (xmax < grid.x0 ? grid.x0 : xmax) 
    return xmin, xmax, ymin, ymax
end

#=
    check_cell_bounds(..., ::NonPeriodicBoundary, ::NonPeriodicBoundary)

Trim cell bounds in the east-west and north-south direction if they exend past
grid due to non-periodic boundary pairs. 
=#
function check_cell_bounds(
    xmin,
    xmax,
    ymin,
    ymax,
    grid,
    ::NonPeriodicBoundary,
    ::NonPeriodicBoundary,
)
    xmin = xmin < grid.x0 ?
        grid.x0 :
        (xmin > grid.xf ? grid.xf : xmin)

    xmax = xmax > grid.xf ?
        grid.xf :
        (xmax < grid.x0 ? grid.x0 : xmax)

    ymin = ymin < grid.y0 ?
        grid.y0 :
        (ymin > grid.yf ? grid.yf : ymin)

    ymax = ymax > grid.yf ?
        grid.yf :
        (ymax < grid.y0 ? grid.y0 : ymax)
    return xmin, xmax, ymin, ymax
end

function center_cell_bounds(
    xidx::Int,
    yidx::Int,
    ns_bound,
    ew_bound,
) where FT
    xmin = (xidx - 1.5) * grid.Δx + grid.x0
    xmax = xmin + grid.Δx
    ymin = (yidx - 1.5) * grid.Δy + grid.y0
    ymax = ymin + grid.Δy
    # Check if cell extends beyond boundaries and if non-periodic, trim it to
    # fit within grid
    return xmin, xmax, ymin, ymax = check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ns_bound,
        ew_bound,
    )
end

#=
    center_cell_poly(...)

Find the coordinates of a given grid cell, centered on a grid line with row yidx
and column xidx. This is offset from the cells within the regular rectilinear
grid by half of a grid cell. Create a polygon that is the shape of this cell.
Inputs:
    xidx        <Int> x index of grid line within list of gridlines (cell column)
    yidx        <Int> y index of grid line within list of gridlines (cell row)
    grid        <RegRectilinearGrid> model's grid 
    ns_bound    <AbstractBoundary> type of either north or south boundary - for
                    checking if periodic
    ew_bound    <AbstractBoundary> type of either east or west boundary - for
                    checking if perioidic

Note that cell bounds will be adjusted depending on if the bounds are
periodic. Cells cannot extend outside of non-periodic boundaries and thus
will be trimmed at boundaries. Therefore, if indices place cell completely
outside of grid, could return a line at the edge of the boundary. 
=#
function center_cell_poly(
    ::Type{FT},
    xidx::Int,
    yidx::Int,
    grid::RegRectilinearGrid,
    ns_bound,
    ew_bound,
) where FT
    xmin, xmax, ymin, ymax = center_cell_bounds(xidx, yidx, grid, ns_bound, ew_bound)
    return _make_bounding_box_polygon(FT, xmin, xmax, ymin, ymax)
end

# Return index as is given non-periodic boundary pair in either x or y direction.
function shift_cell_idx(idx, nlines, ::NonPeriodicBoundary)
    return idx
end

#=
If index is greater than or equal to the grid lines, shift index to equivalent
grid line on opposite side of grid due to periodic boundary. Similarly if given
index is less than 1, shift index to equivalent grid line on opposite side of
grid due to periodic boundary.
For example, the last grid index, nlines, is equivalent to the 1st grid line.
The nlines+1 grid line is equivalent to the 2nd grid line.
=#
function shift_cell_idx(idx, nlines, ::PeriodicBoundary)
    ncells = nlines - 1
    return idx < 1 ? (idx + ncells) : ncells < idx ? (idx - ncells) : idx
end

#-------------- Ocean and Atmosphere on Ice --------------#
"""
    floe_to_grid_info!(...)

Add force from the ice on ocean to ocean force fields (fx & fy) for each grid
cell and update ocean sea ice area fraction (si_area), representing total area
of sea ice in a given cell. Function is called for each sub-floe point.

## _Positional arguments_
- `floeidx::Int`: index of floe within model's floe array
- `xidx::Int`: grid x index that floe's point is within for grid centered on grid lines
- `yidx::Int`: grid column that floe's point is within for grid centered on grid lines
- `τx_ocn::AbstractFloat`: x-stress caused by ocean on point
- `τy_ocn::AbstractFloat`: y-stress caused by ocean on point
- $GRID_DEF
- $DOMAIN_DEF
- `cell_floes::Matrix{CellFloes}`: matrix of `CellFloes`, one for each grid cell
- `ocean::Ocean`: model's ocean
- `coupling_settings::CouplingSettings`: simulation's coupling settings

## _Returns_
- None. Updates `cell_floes` and `scells`.
"""
function floe_to_grid_info!(
    floeidx,
    xidx,
    yidx,
    τx_ocn::FT,
    τy_ocn::FT,
    grid,
    ns_bound,
    ew_bound,
    ocean,
) where {FT}
    # Determine grid cell point is in and if floe is shifted by periodic bounds
    cx = shift_cell_idx(xidx, grid.Nx + 1, ew_bound)
    cy = shift_cell_idx(yidx, grid.Ny + 1, ns_bound)
    n = grid.floe_loc_n_cells[floeidx]
    for k in 1:n                                   # search floe i's own cells
        if grid.floe_loc_cell_x[floeidx, k] == cx && grid.floe_loc_cell_y[floeidx, k] == cy
            ocean.scells_τx[floeidx, k] -= τx_ocn
            ocean.scells_τy[floeidx, k] -= τy_ocn
            ocean.scells_npoints[floeidx, k] += 1
            return
        end
    end
    n += 1                                         # new cell for this floe
    grid.floe_loc_n_cells[floeidx] = n
    grid.floe_loc_cell_x[floeidx, n] = cx
    grid.floe_loc_cell_y[floeidx, n] = cy
    grid.floe_loc_Δx[floeidx, n] = (cx - xidx) * grid.Δx
    grid.floe_loc_Δy[floeidx, n] = (cy - yidx) * grid.Δy
    ocean.scells_τx[floeidx, n] = -τx_ocn
    ocean.scells_τy[floeidx, n] = -τy_ocn
    ocean.scells_npoints[floeidx, n] = 1
    return
end



#=
GPU version:
Calculates subfloe point's cartesian coordinates and index within the grid and 
records them in pre-allocated input arrays, cart_vals and grid_idx. Returns
the index of the last element in cart_vals and grid_idx that corresponds to
elements from given floe.
=#
function calc_subfloe_values!(
    floes::FixedWidthFloes,
    i,
    grid,
    domain,
    cart_vals,
    grid_idx,
)
    # Translate/rotate monte carlo points to floe location/orientation
    α = floes.α[i]
    j = 0 # index in output array
    for k in 1:floes.subfloe_points[i]  # index of monte carlo points
        px = cos(α)*floes.x_subfloe_points[i,k] -
            sin(α)*floes.y_subfloe_points[i,k]  # at origin
        py = sin(α)*floes.x_subfloe_points[i,k] +
            cos(α)*floes.y_subfloe_points[i,k]  # at origin
        x = px + floes.centroid[i,1]  # at centroid
        y = py + floes.centroid[i,2]  # at centroid

        # If point is in bounds, continue to find rest of values
        if in_bounds(x, y, grid, domain.north, domain.east)
            j += 1  # if added to outputs, move to next index in output array
            cart_vals[i, j, 1] = x
            cart_vals[i, j, 2] = y
            grid_idx[i, j, 1], grid_idx[i, j, 2] = find_center_cell_index(
                cart_vals[i, j, 1],
                cart_vals[i, j, 2],
                grid,
            )
        end
    end
    return j # number of point-per-floe inside booundaries
end





function one_way_coupling(
    floes::FixedWidthFloes{FT},
    i, 
    grid, 
    atmos,
    ocean,
    domain,
    coupling_settings,
    consts,
    cart_vals,
    grid_idx,
) where {FT}

    grid.floe_loc_n_cells[i] = 0 # number of ocean cells crossed by each floe

    npoints = calc_subfloe_values!(
        floes,
        i,
        grid,
        domain,
        cart_vals,
        grid_idx,
    )
    if npoints == 0
        floes.tag[i] = remove
    else
        # Add coriolis stress to total stress - same for every point
        xcoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.v[i]
        ycoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.u[i]
        tot_τx = npoints * xcoriolis
        tot_τy = -npoints * ycoriolis
        tot_τtrq = FT(0)
        tot_hflx_factor = FT(0)
        # Determine total stress per-monte carlo point
        for j in 1:npoints
            # Monte carlo point properties
            xcentered = cart_vals[i, j, 1] - floes.centroid[i, 1]
            ycentered = cart_vals[i, j, 2] - floes.centroid[i, 2]
            θ = atan(ycentered, xcentered)
            rad = sqrt(xcentered^2 + ycentered^2)
            upoint = floes.u[i] - floes.ξ[i] * rad * sin(θ)
            vpoint = floes.v[i] + floes.ξ[i] * rad * cos(θ)
            # Stress at monte carlo point from ocean and atmosphere
            τx_atm, τy_atm = atmosphere_forcing(
                cart_vals[i, j, 1],
                cart_vals[i, j, 2],
                upoint,
                vpoint,
                atmos,
                grid,
                domain,
                consts
            )
            τx_ocn, τy_ocn, τx_p∇, τy_p∇, hflx_factor = ocean_forcing(
                cart_vals[i, j, 1],
                cart_vals[i, j, 2],
                upoint,
                vpoint,
                ocean,
                floes.mass[i]/floes.area[i],
                grid,
                domain,
                consts,
            )
            τx = τx_atm + τx_p∇ + τx_ocn
            τy = τy_atm + τy_p∇ + τy_ocn
            # Torque at monte carlo point
            τtrq = (-τx * sin(θ) + τy * cos(θ)) * rad
            # Add values to total stresses
            tot_τx += τx
            tot_τy += τy
            tot_τtrq += τtrq
            tot_hflx_factor += hflx_factor
            # Save floe info onto the grid
            if coupling_settings.two_way_coupling_on 
                floe_to_grid_info!(
                    i,
                    grid_idx[i, j, 1],
                    grid_idx[i, j, 2],
                    τx_ocn,
                    τy_ocn,
                    grid,
                    domain.north,
                    domain.east,
                    ocean,
                )
            end
        end
        # Average forces on ice floe
        floes.fxOA[i] = tot_τx/npoints * floes.area[i]
        floes.fyOA[i] = tot_τy/npoints * floes.area[i]
        floes.trqOA[i] = tot_τtrq/npoints * floes.area[i]
        floes.hflx_factor[i] = tot_hflx_factor/npoints
    end
end

"""
    calc_one_way_coupling!(...)

Preforms calculations needed for one way coupling by calculating floe's forcings
from ocean and atmosphere as well as the heatflux below a given floe.

## _Positional arguments_
- $FLOES_DEF
- $GRID_DEF
- `atmos::Atoms`: model's atmosphere
- `ocean::Ocean`: model's ocean
- $DOMAIN_DEF
- `coupling_settings::CouplingSettings`: simulation coupling settings
- $CONSTS_DEF

## _Returns_
- None. Update each floe's forces, torque, and heatflux factor from ocean/atmosphere.
"""
function calc_one_way_coupling!(
    dev_floes::FixedWidthFloes{FT},
    grid,
    atmos,
    ocean,
    domain,
    coupling_settings,
    consts,
    backend,
) where {FT}
    nfloes = length(dev_floes.height)
    max_points = size(dev_floes.x_subfloe_points, 2)
    cart_vals = adapt(backend, Array{FT,3}(undef, nfloes, max_points, 2))
    grid_idx  = adapt(backend, Array{Int,3}(undef, nfloes, max_points, 2))  

    launch_per_floe!(one_way_coupling,
        backend,
        dev_floes,
        grid,
        atmos,
        ocean,
        domain,
        coupling_settings,
        consts,
        cart_vals,
        grid_idx,
    )
end



function ice_area_in_cells!(floes::FixedWidthFloes, i, grid, ocean, domain, B, S)
    for k in 1:grid.floe_loc_n_cells[i]
        cx = grid.floe_loc_cell_x[i, k]
        cy = grid.floe_loc_cell_y[i, k]
        xmin, xmax, ymin, ymax = center_cell_bounds(cx, cy, grid, domain.north, domain.east)
        A = _area_in_rectangle(floes, i, grid.floe_loc_Δx[i, k], grid.floe_loc_Δy[i, k],
            xmin, xmax, ymin, ymax)
        ocean.scells_area[i, k] = A
        if A > 0
            c = cx + (cy - 1) * (grid.Nx + 1)
            ocean.scells_key[(i - 1) * S + k] = (c - 1) * B + (i - 1) * S + k
        end
    end
    return
end

function ocean_stress_in_cell!(ocean, ix, iy, grid, atmos, B, S, consts, ρi, Δt)
    FT = eltype(ocean.τx)
    keys = ocean.scells_key
    c = ix + (iy - 1) * (grid.Nx + 1)
    τx, τy, si_frac = zero(FT), zero(FT), zero(FT)
    idx = _first_key_at_least(keys, (c - 1) * B + 1)
    while idx <= length(keys) && keys[idx] <= c * B
        r = keys[idx] - 1
        i, k = (r % B) ÷ S + 1, r % S + 1
        A = ocean.scells_area[i, k]
        τx += (ocean.scells_τx[i, k] / ocean.scells_npoints[i, k]) * A
        τy += (ocean.scells_τy[i, k] / ocean.scells_npoints[i, k]) * A
        si_frac += A
        idx += 1
    end
    if si_frac > 0
        τx /= si_frac
        τy /= si_frac
        si_frac /= grid.Δx * grid.Δy
    end
    Δu = atmos.u[ix, iy] - ocean.u[ix, iy]
    Δv = atmos.v[ix, iy] - ocean.v[ix, iy]
    norm = sqrt(Δu^2 + Δv^2)
    ocean.τx[ix, iy] = τx + consts.ρa * consts.Cd_ao * (1 - si_frac) * norm * Δu
    ocean.τy[ix, iy] = τy + consts.ρa * consts.Cd_ao * (1 - si_frac) * norm * Δv
    ocean.si_frac[ix, iy] = si_frac
    ocean.hflx_factor[ix, iy] = Δt * consts.k / (ρi * consts.L) *
        (ocean.temp[ix, iy] - atmos.temp[ix, iy])
    return
end

# first index of sorted `keys` with keys[idx] ≥ v (length(keys) + 1 if none)
function _first_key_at_least(keys, v)
    lo, hi = 1, length(keys) + 1
    while lo < hi
        mid = (lo + hi) >>> 1
        keys[mid] < v ? (lo = mid + 1) : (hi = mid)
    end
    return lo
end

"""
    calc_two_way_coupling!(...)

Calculate effects of ice and atmosphere on the ocean and update ocean stress
fields and sea ice fraction.

## _Positional arguments_
- $FLOES_DEF
- $GRID_DEF
- `atmos::Atoms`: model's atmosphere
- `ocean::Ocean`: model's ocean
- $DOMAIN_DEF
- $FLOE_SETTINGS_DEF
- $CONSTS_DEF
- $ΔT_DEF

## _Returns_
- None. Update's ocean's stress fields and heatflux factor field. 
"""
function calc_two_way_coupling!(
    dev_floes::FixedWidthFloes{FT},
    i,
    grid,
    atmos,
    ocean,
    domain,
    floe_settings,
    consts,
    Δt,
    backend,
) where {FT}

    launch_per_floe!(two_way_coupling,
        backend,
        dev_floes,
        grid,
        atmos,
        ocean,
        domain,
        floe_settings,
        consts,
        Δt,
        backend,
    )
    return
end

function calc_two_way_coupling!(
    dev_floes::FixedWidthFloes,
    grid,
    atmos,
    ocean,
    domain,
    floe_settings,
    consts,
    Δt,
    backend
)
    S = size(ocean.scells_area, 2)
    B = length(ocean.scells_key)
    fill!(ocean.scells_key, typemax(Int))                                  # stage 0
    launch_per_floe!(ice_area_in_cells!, backend, dev_floes, grid, ocean, domain, B, S)  # stage 1
    sort_keys!(ocean.scells_key, backend)                                   # stage 2 (section 6)
    launch_per_cell!(ocean_stress_in_cell!, backend, ocean, grid, atmos, B, S,
        consts, floe_settings.ρi, Δt)                                       # stage 3
    return
end

"""
    timestep_coupling!(...)

Calculates the effects of the ocean and atmosphere on the ice and the effects of
the ice and atmosphere on the ocean if the coupling is two-way.

## _Positional arguments_
- $MODEL_DEF
- $ΔT_DEF
- $CONSTS_DEF
- `coupling_settings::CouplingSettings`:: simulation coupling settings
- $FLOE_SETTINGS_DEF

## _Returns_
- None. Updates each floe's ocean/atmosphere forcings (fxOA, fyOA, torqueOA)
    and calculates stresses on each ocean grid cell from ice and atmosphere if
    two-way coupling is enabled in coupling_settings       
"""
function timestep_coupling!(
    dev_floes,
    dev_model,
    model,
    Δt,
    consts,
    coupling_settings,
    floe_settings,
    backend,
)
    empty!.(model.grid.floe_locations)
    if coupling_settings.two_way_coupling_on
        empty!.(model.ocean.scells)
    end
    calc_one_way_coupling!(
        dev_floes,
        dev_model.grid,
        dev_model.atmos,
        dev_model.ocean,
        dev_model.domain,
        coupling_settings,
        consts,
        backend,
    )
    if coupling_settings.two_way_coupling_on
        calc_two_way_coupling!(
            dev_floes,
            dev_model.grid,
            dev_model.atmos,
            dev_model.ocean,
            dev_model.domain,
            floe_settings,
            consts,
            Δt,
            backend,
        )
    end
    return
end
