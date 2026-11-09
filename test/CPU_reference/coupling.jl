using Subzero.Interpolations: linear_interpolation
#=
Reference implementation of `timestep_coupling!`, copied from the CPU version on
`main` so that GPU/kernel versions can be checked against it.
=#

# lazily access floe within struct array floe field
_ref_get_floe(floes::StructArray, i::Int) = LazyRow(floes, i)

# Return index as is given non-periodic boundary pair in either x or y direction.
function _ref_shift_cell_idx(idx, nlines, ::Subzero.NonPeriodicBoundary)
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
function _ref_shift_cell_idx(idx, nlines, ::PeriodicBoundary)
    ncells = nlines - 1
    return idx < 1 ? (idx + ncells) : ncells < idx ? (idx - ncells) : idx
end

#=
With all non-periodic boundaries, points outside of the grid in both the x and y
are defined to be out of bounds since these points can't be interpolated as we
don't have any information on the ocean outside of the grid.

Returns a Boolean that is true if both xr and yr are within domain boundaries 
and false otherwise.
=#
function _ref_in_bounds(
    xr,
    yr,
    grid,
    ::Subzero.NonPeriodicBoundary,
    ::Subzero.NonPeriodicBoundary,
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
function _ref_in_bounds(
    xr,
    yr,
    grid,
    ::Subzero.NonPeriodicBoundary,
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
function _ref_in_bounds(
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::Subzero.NonPeriodicBoundary,
)
    return grid.x0 <= xr <= grid.xf
end

#=
With all periodic boundaries, all points are considered to be in-bounds.
Returns a Boolean that is true regardless of point values.
=#
function _ref_in_bounds(
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::PeriodicBoundary,
)
    return true
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
function _ref_find_Subzero.center_cell_index(xp, yp, grid::Union{RegRectilinearGrid, dev_RegRectilinearGrid})
    xidx = floor(Int, (xp - grid.x0)/(grid.Δx) + 0.5) + 1
    yidx = floor(Int, (yp - grid.y0)/(grid.Δy) + 0.5) + 1
    return xidx, yidx
end

#=
Add floe information to grid cell objects (cfloes and scell) that floe sits within. 

This function dispatch is called when two-way coupling is on. Thus, in addition to recording
which floes are within which grid cell, it also records the stresses from the ice onto the ocean
within each grid cell. 

Floe information is saved within the Subzero.CellFloes (cfloes) object, which keeps track of a list of floes within a
given grid cell. If the booundaries are periodic and a floe's centroid is on the opposite side of
of the domain part of its shape, then a Δx and Δy are recorded to note that offset. 

The Subzero.CellStresses (scell) object aggragates the stresses from ice on ocean within a grid cell
from each floes' sub-floe points. τx and τy are the x-directional and y-directional stress
from sub-floe points on the ocean. 
=#
function _ref_add_point!(
    cfloes::Subzero.CellFloes,
    scell::Subzero.CellStresses,
    floeidx,
    τx,
    τy,
    Δx,
    Δy,
)
    if isempty(cfloes.floeidx) || cfloes.floeidx[end] != floeidx
        push!(cfloes.floeidx, floeidx)
        push!(cfloes.Δx, Δx)
        push!(cfloes.Δy, Δy)
        push!(scell.τx, τx)
        push!(scell.τy, τy)
        push!(scell.npoints, 1)
    else
        scell.τx[end] += τx
        scell.τy[end] += τy
        scell.npoints[end] += 1
    end
    return
end

#=
Add floe information to grid cell objects (cfloes and scell) that floe sits within. 

This function dispatch is called when one-way coupling is on. Thus, informatin on which floes are in which
grid cells is recorded. 

Floe information is saved within the Subzero.CellFloes (cfloes) object, which keeps track of a list of floes within a
given grid cell. If the booundaries are periodic and a floe's centroid is on the opposite side of
of the domain part of its shape, then a Δx and Δy are recorded to note that offset. 
=#
function _ref_add_point!(
    cfloes::Subzero.CellFloes,
    floeidx,
    Δx,
    Δy,
)
    if isempty(cfloes.floeidx) || cfloes.floeidx[end] != floeidx
        push!(cfloes.floeidx, floeidx)
        push!(cfloes.Δx, Δx)
        push!(cfloes.Δy, Δy)
    end
    return
end



#=
CPU function that calculates subfloe point's cartesian coordiantes and index
within the grid and records them in pre-allocated input arrays, cart_vals 
and grid_idx. Returns the index of the last element in cart_vals and grid_idx
that corresponds to elements from given floe.
=#
function _ref_calc_subfloe_values!(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    grid,
    domain,
    cart_vals,
    grid_idx,
) where {FT<:AbstractFloat}
    # Translate/rotate monte carlo points to floe location/orientation
    α = floe.α
    j = 0  # index in output array
    for i in eachindex(floe.x_subfloe_points)  # index of monte carlo points
        px = cos(α)*floe.x_subfloe_points[i] -
            sin(α)*floe.y_subfloe_points[i]  # at origin
        py = sin(α)*floe.x_subfloe_points[i] +
            cos(α)*floe.y_subfloe_points[i]  # at origin
        x = px + floe.centroid[1]  # at centroid
        y = py + floe.centroid[2]  # at centroid
        # If point is in bounds, continue to find rest of values
        if _ref_in_bounds(x, y, grid, domain.north, domain.east)
            j += 1  # if added to outputs, move to next index in output array
            cart_vals[j, 1] = x
            cart_vals[j, 2] = y
            grid_idx[j, 1], grid_idx[j, 2] = _ref_find_Subzero.center_cell_index(
                cart_vals[j, 1],
                cart_vals[j, 2],
                grid,
            )
        end
    end
    return j  # last spot in output arrays that has values for this floe
end

#-------------- Interpolation of Ocean and Atmosphere --------------#

#=
    find_interp_knots(..., find_interp_knots)

Find indicies in list of grid lines that surround points with indicies
'point_idx', with a buffer of Δd indices on each side of the points. In this
case, the points are being considered near a periodic boundary, which means that
they can loop around to the other side of the grid. If these points exist, we
extend the grid lines to cover the points and buffer. 
Inputs:
    point_idx   <Vector{Int}> vector of point indices representing the grid line
                    they are nearest
    ncells      <Int> number of grid cells in given dimension
    glines      <Vector or Range> grid line values
    L           <AbstractFloat> length of grid in given dimension
    Δd          <Int> number of buffer grid cells to include on either side of
                    the provided indicies 
                <PeriodicBoundary> dispatching on periodic boundary
Outputs:
    knots       <Vector{AbstractFloat}> interpolation knots - grid line values
    knot_idx    <Vector{Int}> - indices of grid line values within grid list of
                    grid lines.
Note:
    The grid values are extended if points expand past gridlines, however, the
    indices are within the grid. For example, consider a grid where the maximum
    grid value is 1e5, with grid cells of length 1e4. One of the knot values
    might be 1.1e5, however, its index would be 2 since the grid line at 1e5 is
    equivalent to the first grid line since it is periodic, and 1.1e5 is one
    grid cell length past that value.

    This function depends on the ocean being periodic in the given direction.
    We assume that first grid line and the last grid line are the same, and have
    the same values within the ocean/atmosphere. These are not repeated in the
    knots, but rather only one is used. So if there are 10 grid lines, grid line
    1 and 10 are the equivalent and we use grid line 1 exclusively. 
=#
function _ref_find_interp_knots(
    point_idx,
    ncells,
    glines,
    L,
    Δd::Int,
    ::PeriodicBoundary,
)
    min_line, max_line = extrema(point_idx)
    # Grid lines surrounding points with buffers
    # Point close to ith grid line could be between the i and i-1 grid line
    min_line -= (Δd + 1)
    # Point close to ith grid line could be between the i and i+1 grid line
    max_line += (Δd + 1)

    # Find out-of-bounds (oob) indices and the in-bounds (within grid) indices
    # Out of bounds on south or west side of domain 
    low_oob_idx = 1:0  # empty range
    # Out of bounds on north or east side of domain 
    high_oob_idx = 1:0  # empty range
    in_bounds_idx = 
        # Last gird line is equal to first grid line
        if min_line < 1 && max_line > ncells
            low_oob_idx = (min_line + ncells):ncells
            high_oob_idx = 1:(max_line-ncells)
            1:ncells
        elseif min_line < 1
            low_oob_idx = (min_line + ncells):ncells
            1:max_line
        elseif max_line > ncells
            high_oob_idx = 1:(max_line-ncells)
            min_line:ncells
        else
            min_line:max_line
        end
    knot_idx = [low_oob_idx; in_bounds_idx; high_oob_idx]
    #= Adjust out-of-bound values by grid length so there isn't a jump in
    interpolation spacing =#
    knots = [
        glines[low_oob_idx] .- L;
        glines[in_bounds_idx];
        glines[high_oob_idx] .+ L
    ]

    return knots, knot_idx
end

#=
    find_interp_knots(..., ::Subzero.NonPeriodicBoundary)

Find indicies in list of grid lines that surround points with indicies
'point_idx' with a buffer of Δd indices on each side of the points. In this
case, the points are being considered near a NON-periodic boundary, so we cut
off the possible indices past the edge of the grid. 
Inputs:
    point_idx   <Vector{Int}> vector of indices representing the grid line they
                    are nearest
    ncells      <Int> number of grid cells in given dimension
    glines      <Vector or Range> grid line values
    L           <AbstractFloat> length of grid in given dimension
    Δd          <Int> number of buffer grid cells to include on either side of
                    the provided indicies 
                <PeriodicBoundary> dispatching on periodic boundary
Outputs:
    knots       <Vector{AbstractFloat}> interpolation knots - grid line values
    knot_idx    <Vector{Int}> - indices of grid line values within grid list of
                    grid lines.
Note:
    Only knots within the grid will be returned since this is a non-periodic
    boundary.
=#
function _ref_find_interp_knots(
    point_idx,
    ncells,
    glines,
    L,
    Δd::Int,
    ::Subzero.NonPeriodicBoundary,
)
    nlines = ncells + 1
    min_line, max_line = extrema(point_idx)
    # point close to ith grid line could be between the i and i-1 grid line
    min_line -= (Δd + 1)
    # point close to ith grid line could be between the i and i+1 grid line
    max_line += (Δd + 1)
    # Can't be outside grid given unknown ocean/atmosphere conditions
    min_line = (min_line < 1) ? 1 : min_line
    max_line = (max_line > nlines) ? nlines : max_line
    return glines[min_line:max_line], [min_line:max_line;]
end

#=
    _ref_mc_interpolation(
        npoints,
        mc_grid_idx,
        grid,
        domain,
        atmos,
        ocean,
        coupling_settings,
    )

Create and returns interpolation objects for atmosphere u and v velocities, and
ocean u and v velocities, in addition to ocean's heatflux factor.
Inputs:
    npoints             <Int> number of monte carlo points to consider - the
                            number of rows to use in mc_cart and mc_grid_idx
    grid_idx            <Matrix{Int}> index of monte carlo points within the
                            grid - nx2 matrix of indices where the first column
                            is the grid column index and the second column is
                            the grid row index for cells centered on grid lines
    grid                <AbstractRectilinearGrid> model grid
    domain              <Domain> model domain
    atmos               <Atmos> model atmosphere
    ocean               <Ocean> model ocean
    coupling_settings   <CouplingSettings> simulation's coupling settings
Outputs:
    uatm_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the atompshere u velocity onto point
    vatm_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the atompshere v velocity onto point
    uocn_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean u velocity onto point
    vocn_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean v velocity onto point
    hflx_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean heatflux factor velocity onto point
=#
function _ref_mc_interpolation(
    npoints,
    grid_idx,
    grid,
    domain,
    atmos,
    ocean,
    coupling_settings,
)
    xidx = @view grid_idx[1:npoints, 1]
    yidx = @view grid_idx[1:npoints, 2]

    # Find knots and indices of knots for monte carlo interpolation
    xknots, xknot_idx = _ref_find_interp_knots(
        xidx,
        grid.Nx,
        grid.x0:grid.Δx:grid.xf,
        grid.xf - grid.x0,
        coupling_settings.Δd,
        domain.east,
    )
    yknots, yknot_idx = _ref_find_interp_knots(
        yidx,
        grid.Ny,
        grid.y0:grid.Δy:grid.yf,
        grid.yf - grid.y0,
        coupling_settings.Δd,
        domain.north,
    )

    knots = (xknots, yknots)

    # Atmos Interpolation objects for Monte Carlo Points
    uatm_interp = linear_interpolation(
        knots,
        @view(atmos.u[xknot_idx, yknot_idx]),
    )
    vatm_interp = linear_interpolation(
        knots,
        @view(atmos.v[xknot_idx, yknot_idx]),
    )
    
    # Ocean Interpolation objects for Monte Carlo Points
    uocn_interp = linear_interpolation(
        knots,
        @view(ocean.u[xknot_idx, yknot_idx]),
    )
    vocn_interp = linear_interpolation(
        knots,
        @view(ocean.v[xknot_idx, yknot_idx]),
    )
    hflx_interp = linear_interpolation(
        knots,
        @view(ocean.hflx_factor[xknot_idx, yknot_idx]),
    )

    return uatm_interp, vatm_interp, uocn_interp, vocn_interp, hflx_interp
end



#-------------- Ocean and Atmosphere on Ice --------------#
#=
    calc_atmosphere_forcing(...)

Calculates the stresses on a floe from the atmosphere above at given monte
carlo point.

## _Positional arguments_
- `xr::AbstractFloat`: x-coordiantes of points to interplate on
- `yr::AbstractFloat`: y-coordiantes of points to interplate on
- `upoint::AbstractFloat`: u velocity of floe at point
- `vpoint::AbstractFloat`: v velocity of floe at point
- `uatm_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the atompshere u velocity onto point
- `vatm_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the atompshere v velocity onto point
- $CONSTS_DEF

## _Returns_
- `τx_atm::AbstractFloat`: stress from atmosphere on floe in x-direction at given point
- `τy_atm::AbstractFloat`: stress from atmosphere on floe iny-direction at given point
=# 
function _ref_calc_atmosphere_forcing(
    xr, 
    yr,
    upoint,
    vpoint,
    uatm_interp,
    vatm_interp,
    c,
)
    # Atmosphere velocities at monte carlo point
    uatm = uatm_interp(xr, yr) 
    vatm = vatm_interp(xr, yr) 

    # Stress on ice from atmopshere
    Δu_AI = uatm - upoint
    Δv_AI = vatm - vpoint
    norm = sqrt(Δu_AI^2 + Δv_AI^2)
    τx_atm = c.ρa * c.Cd_ia * norm * Δu_AI
    τy_atm = c.ρa * c.Cd_ia * norm * Δv_AI
    return τx_atm, τy_atm
end


#=
    calc_ocean_forcing!(...)

Calculates the stresses on a floe from the ocean above at given monte carlo
point.

## _Positional arguments_
- `xr::AbstractFloat`: x-coordiantes of points to interplate on
- `yr::AbstractFloat`: y-coordiantes of points to interplate on
- `upoint::AbstractFloat`: u velocity of floe at point
- `vpoint::AbstractFloat`: v velocity of floe at point
- `uocn_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the ocean u velocity onto point
- `vocn_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the ocean v velocity onto point
- `hflx_interp::InterplationsObject`: linear interpolation function from
    Interpolations.jl that takes in two arguments (x, y) and interpolates
    the ocean heatflux factor velocity onto point
- `ma_ratio::AbstractFloat`: floe's mass to area ratio
- $CONSTS_DEF

## _Returns_
- `τx_ocn::AbstractFloat`: stress from ocean velocity on floe in x-direction at given point
- `τy_ocn::AbstractFloat`: stress from ocean velocity on floe in y-direction at given point
- `τx_pressure∇::AbstractFloat`: stress from ocean pressure gradient on floe in x-direction at given point
- `τy_pressure∇::AbstractFloat`: stress from ocean pressure gradient on floe in y-direction at given point
- `hflx_factor::AbstractFloat`: heatflux factor at given point from the heatflux factors of ocean below floe
=#
function _ref_calc_ocean_forcing!(
    xr,
    yr,
    upoint,
    vpoint,
    uocn_interp,
    vocn_interp,
    hflx_interp,
    ma_ratio,
    c,  # constants
)
    uocn = uocn_interp(xr, yr)
    vocn = vocn_interp(xr, yr)
    hflx_factor = hflx_interp(xr, yr)
    Δu_OI = uocn - upoint
    Δv_OI = vocn - vpoint
    norm = sqrt(Δu_OI^2 + Δv_OI^2)
    τx_ocn = c.ρo*c.Cd_io * norm * (cos(c.turnθ) * Δu_OI - sin(c.turnθ) * Δv_OI)
    τy_ocn = c.ρo*c.Cd_io * norm * (sin(c.turnθ) * Δu_OI + cos(c.turnθ) * Δv_OI)
    τx_pressure∇ = -ma_ratio * c.f * vocn
    τy_pressure∇ = ma_ratio * c.f * uocn
    return τx_ocn, τy_ocn, τx_pressure∇, τy_pressure∇, hflx_factor
end


#=
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
- `cell_floes::Matrix{Subzero.CellFloes}`: matrix of `Subzero.CellFloes`, one for each grid cell
- `scells::Matrix{Subzero.CellStresses}`: matrix of `Subzero.CellStressess`, one for each grid cell
- `coupling_settings::CouplingSettings`: simulation's coupling settings

## _Returns_
- None. Updates `cell_floes` and `scells`.
=#
function _ref_floe_to_grid_info!(
    floeidx,
    xidx,
    yidx,
    τx_ocn::FT,
    τy_ocn::FT,
    grid,
    ns_bound,
    ew_bound,
    scells,
    coupling_settings,
) where {FT}
    # Determine grid cell point is in and if floe is shifted by periodic bounds
    shifted_xidx = _ref_shift_cell_idx(xidx, grid.Nx + 1, ew_bound)
    shifted_yidx = _ref_shift_cell_idx(yidx, grid.Ny + 1, ns_bound)
    Δx = (shifted_xidx - xidx) * (grid.Δx)
    Δy = (shifted_yidx - yidx) * (grid.Δy)
    if coupling_settings.two_way_coupling_on 
        # If two-way coupling, save stress on ocean per cell
        _ref_add_point!(
            grid.floe_locations[shifted_xidx, shifted_yidx],
            scells[shifted_xidx, shifted_yidx],
            floeidx,
            -τx_ocn,
            -τy_ocn,
            Δx,
            Δy,
        )
    else
        _ref_add_point!(
            grid.floe_locations[shifted_xidx, shifted_yidx],
            floeidx,
            Δx,
            Δy,
        )
    end
    return
end



#=
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
=#
function _ref_calc_one_way_coupling!(
    floes::StructArray{Floe{FT}},
    grid,
    atmos,
    ocean,
    domain,
    coupling_settings,
    consts,
) where {FT}
    max_points = maximum(length, floes.x_subfloe_points)
    cart_vals = Matrix{FT}(undef, max_points, 2)
    grid_idx = Matrix{Int}(undef, max_points, 2)
    for i in eachindex(floes)
        # Monte carlo point cartesian coordinates and grid cell indices
        npoints = _ref_calc_subfloe_values!(
            _ref_get_floe(floes, i),
            grid,
            domain,
            cart_vals,
            grid_idx,
        )
        if npoints == 0
            floes.status[i].tag = Subzero.remove
        else
            # Interpolaters for ocean and atmosphere
            uatm_int, vatm_int, uocn_int, vocn_int, hflx_int = _ref_mc_interpolation(
                npoints,
                grid_idx,
                grid,
                domain,
                atmos,
                ocean,
                coupling_settings,
            )

            # Add coriolis stress to total stress - same for every point
            xcoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.v[i]
            ycoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.u[i]
            tot_τx = npoints * xcoriolis
            tot_τy = -npoints * ycoriolis
            tot_τtrq = FT(0)
            tot_hflx_factor = FT(0)
            ma_ratio = floes.mass[i]/floes.area[i]
            # Determine total stress per-monte carlo point
            for j in 1:npoints
                # Monte carlo point properties
                xcentered = cart_vals[j, 1] - floes.centroid[i][1]
                ycentered = cart_vals[j, 2] - floes.centroid[i][2]
                θ = atan(ycentered, xcentered)
                rad = sqrt(xcentered^2 + ycentered^2)
                upoint = floes.u[i] - floes.ξ[i] * rad * sin(θ)
                vpoint = floes.v[i] + floes.ξ[i] * rad * cos(θ)
                # Stress at monte carlo point from ocean and atmosphere
                τx_atm, τy_atm = _ref_calc_atmosphere_forcing(
                    cart_vals[j, 1], 
                    cart_vals[j, 2],
                    upoint,
                    vpoint,
                    uatm_int,
                    vatm_int,
                    consts,
                )
                τx_ocn, τy_ocn, τx_p∇, τy_p∇, hflx_factor = _ref_calc_ocean_forcing!(
                    cart_vals[j, 1],
                    cart_vals[j, 2],
                    upoint,
                    vpoint,
                    uocn_int,
                    vocn_int,
                    hflx_int,
                    ma_ratio,
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
                _ref_floe_to_grid_info!(
                    i,
                    grid_idx[j, 1],
                    grid_idx[j, 2],
                    τx_ocn,
                    τy_ocn,
                    grid,
                    domain.north,
                    domain.east,
                    ocean.scells,
                    coupling_settings,
                )
            end
            # Average forces on ice floe
            floes.fxOA[i] = tot_τx/npoints * floes.area[i]
            floes.fyOA[i] = tot_τy/npoints * floes.area[i]
            floes.trqOA[i] = tot_τtrq/npoints * floes.area[i]
            floes.hflx_factor[i] = tot_hflx_factor/npoints
        end
    end
end


#=
    _ref_calc_two_way_coupling!(...)

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
=#
function _ref_calc_two_way_coupling!(
    floes::StructArray{Floe{FT}},
    grid::RegRectilinearGrid,
    atmos,
    ocean,
    domain,
    floe_settings,
    consts,
    Δt,
) where {FT}
    # Determine force from floe on each grid cell it is in
    cell_area = grid.Δx * grid.Δy
    Threads.@threads for cartidx in CartesianIndices(ocean.scells)
        ocean.τx[cartidx] = FT(0)
        ocean.τy[cartidx] = FT(0)
        ocean.si_frac[cartidx] = FT(0)
        τocn = ocean.scells[cartidx]
        floe_locations = grid.floe_locations[cartidx]
        if !isempty(floe_locations.floeidx)
            # Coordinates of grid cell
            cell_poly = Subzero.center_cell_poly(
                FT,
                cartidx[1],
                cartidx[2],
                grid,
                domain.north,
                domain.east
            )
            for i in eachindex(floe_locations.floeidx)
                floe_poly = Subzero._translate_poly(FT,
                    floes.poly[floe_locations.floeidx[i]],
                    floe_locations.Δx[i],
                    floe_locations.Δy[i],
                )::Subzero.Polys{FT}
                floe_area_in_cell = sum(
                    Subzero.area_poly.(Subzero.intersect_polys(cell_poly, floe_poly), FT)
                )
                if floe_area_in_cell > 0
                    # Add forces and area to ocean fields
                    ocean.τx[cartidx] += (τocn.τx[i]/τocn.npoints[i]) * floe_area_in_cell
                    ocean.τy[cartidx] += (τocn.τy[i]/τocn.npoints[i]) * floe_area_in_cell
                    ocean.si_frac[cartidx] += floe_area_in_cell
                end
            end
            if ocean.si_frac[cartidx] > 0
                # Divide by total floe area in cell to get ocean stress
                ocean.τx[cartidx] /= ocean.si_frac[cartidx]
                ocean.τy[cartidx] /= ocean.si_frac[cartidx]
                # Divide by cell area to get sea-ice fraction
                ocean.si_frac[cartidx] /= cell_area
            end
        end
        Δu_AO = atmos.u[cartidx] - ocean.u[cartidx]
        Δv_AO = atmos.v[cartidx] - ocean.v[cartidx]
        ocn_frac = 1 - ocean.si_frac[cartidx]
        norm = sqrt(Δu_AO^2 + Δv_AO^2)
        ocean.τx[cartidx] += consts.ρa * consts.Cd_ao * ocn_frac * norm * Δu_AO
        ocean.τy[cartidx] += consts.ρa * consts.Cd_ao * ocn_frac * norm * Δv_AO
        # Not sure this is where the heatflux should be??
        ocean.hflx_factor[cartidx] = Δt * consts.k/(floe_settings.ρi*consts.L) *
            (ocean.temp[cartidx] - atmos.temp[cartidx])
    end
    return
end

#=
    _ref_timestep_coupling!(...)

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
=#
function _ref_timestep_coupling!(
    model,
    Δt,
    consts,
    coupling_settings,
    floe_settings,
)
    empty!.(model.grid.floe_locations)
    if coupling_settings.two_way_coupling_on
        empty!.(model.ocean.scells)
    end
    _ref_calc_one_way_coupling!(
        model.floes,
        model.grid,
        model.atmos,
        model.ocean,
        model.domain,
        coupling_settings,
        consts,
    )
    if coupling_settings.two_way_coupling_on
        _ref_calc_two_way_coupling!(
            model.floes,
            model.grid,
            model.atmos,
            model.ocean,
            model.domain,
            floe_settings,
            consts,
            Δt,
        )
    end
    return
end