"""
Functions needed for coupling the ice, ocean, and atmosphere.
"""

"""
    AbstractSubFloePointsGenerator

Abstract type for parameters determining generation of sub-floe points used for
interpolation. The points generated using these parameters will be used to find
stresses on each floe from the ocean and the atmosphere. There must be a
`generate_subfloe_points` function that dispatches off of the subtype of
AbstractSubFloePointsGenerator to generate the points for a given floe.
Points generated must all be within a given floe.
"""
abstract type AbstractSubFloePointsGenerator{FT<:AbstractFloat} end

"""
    MonteCarloPointsGenerator

Subtype of AbstractSubFloePointsGenerator that defines parameters needed to
generate a set of random monte carlo points within a given floe. `npoints` is
the number of monte carlo points to generate within the floe's bounding box -
the floe will not end up with this many points as all points outside of the floe
will be removed. `ntries` is the number of tries to generate a set of points
within the floe that have a smaller error than `err`.
"""
@kwdef struct MonteCarloPointsGenerator{
    FT <: AbstractFloat,
} <: AbstractSubFloePointsGenerator{FT}
    npoints::Int = 1000
    ntries::Int = 10
    err::FT = 0.1

    function MonteCarloPointsGenerator{FT}(
        npoints,
        ntries,
        err,
    ) where {FT <: AbstractFloat}
        if npoints < 1
            throw(ArgumentError("Interpolation cannot be preformed with no \
            monte carlo points. Field npoints must be positive."))
        end
        if ntries < 1
            throw(ArgumentError("Monte carlo points cannot be generated without\
             trying at least once. Field ntried must be positive."))
        end
        if err < 0 || err > 1
            throw(ArgumentError("Field err must be between 0 and 1."))
        end
        return new{FT}(npoints, ntries, err)
    end
end

"""
    MonteCarloPointsGenerator(::Type{FT}; kwargs...)

A float type FT can be provided as the first argument of any
MonteCarloPointsGenerator constructor. A MonteCarloPointsGenerato of type FT
will be created by passing all other arguments to the correct constructor. 
"""
MonteCarloPointsGenerator(
    ::Type{FT},
    args...;
    kwargs...,
) where {FT <: AbstractFloat} =
    MonteCarloPointsGenerator{FT}(args...; kwargs...)

"""
    MonteCarloPointsGenerator(; kwargs...)

If type isn't specified, MonteCarloPointsGenerator(; kwargs...) will be of type
Float64 and the correct constructor will be called with all other arguments.
"""
MonteCarloPointsGenerator(args...; kwargs...) =
    MonteCarloPointsGenerator{Float64}(args...; kwargs...)

"""
    SubGridPointsGenerator

Subtype of AbstractSubFloePointsGenerator that defines parameters needed to
generate a set of points on a "subgrid" within the floe where the subgrid is a
regular rectilinar grid with cells of size `Δg` in both width and height. If
two-way coupling, `Δg` should be smaller than the grid's Δx and Δy so that there
is at least one point in each grid cell that the floe occupies.
"""
@kwdef struct SubGridPointsGenerator{
    FT <: AbstractFloat,
} <: AbstractSubFloePointsGenerator{FT}
    Δg::FT

    function SubGridPointsGenerator{FT}(Δg) where {FT <: AbstractFloat}
        if Δg <= 0
            throw(ArgumentError("Field Δg must be positive as it is the width \
            and height value of the sub-floe grid cells."))
        end
        return new{FT}(Δg)
    end
end

"""
    SubGridPointsGenerator(::Type{FT}; kwargs...)

A float type FT can be provided as the first argument of any
SubGridPointsGenerator constructor. A SubGridPointsGenerator of type FT will be
created by passing all other arguments to the correct constructor. 
"""
SubGridPointsGenerator(
    ::Type{FT},
    args...;
    kwargs...,
) where {FT <: AbstractFloat} =
    SubGridPointsGenerator{FT}(args...; kwargs...)

"""
    SubGridPointsGenerator(; kwargs...)

If type isn't specified, SubGridPointsGenerator(; kwargs...) will be of type
Float64 and the correct constructor will be called with all other arguments.
"""
SubGridPointsGenerator(args...; kwargs...) =
    SubGridPointsGenerator{Float64}(args...; kwargs...)

"""
    SubGridPointsGenerator{FT}(grid, npoint_per_cell)

SubGridPointsGenerator constructor that uses the simulation grid and the desired
number of sub-floe points per simulation grid cell to determine the correct
value of Δg
Inputs:
    grid                <RegRectilinearGrid> simulation's grid
    npoint_per_cell     <Int> number of points per grid cell, where the grid
                            is redefined to have width and height equal to the
                            minimum of Δx and Δy
Output:
    SubGridPointsGenerator with Δg defined to be the minimum of Δx and Δy over
    the number of desired points per grid cell
"""
SubGridPointsGenerator{FT}(
    grid::RegRectilinearGrid,
    npoint_per_cell::Int,
) where {FT <: AbstractFloat} =
    SubGridPointsGenerator{FT}(
        min(grid.Δx, grid.Δy) / npoint_per_cell / sqrt(2)
    )

"""
    generate_subfloe_points(
        point_generator
        poly,
        centroid,
        area,
        status,
        rng
    )

Generate monte carlo points centered on the origin within the floe according to
parameters defined in the point_generator argument.
Inputs:
    point_generator     <MonteCarloPointsGenerator> monte carlo point generator
    poly                <Polygon> Polygon representing floe shape
    centroid            <Matrix> floe's centroid
    area                <AbstractFloat> floe's area
    status              <Status> floe status (i.e. active, fuse in simulation)
    rng                 <AbstractRNG> random number generator to generate monte
                            carlo points
Ouputs:
    x_sub_floe  <Vector{FT}> vector of sub-floe grid points x-coords within floe
    y_sub_floe  <Vector{FT}> vector of sub-floe grid points y-coords within floe
    status      <Status> floe's status post generation, changed to remove if 
                    generation is unsuccessful
"""
function generate_subfloe_points(
    point_generator::MonteCarloPointsGenerator{FT},
    poly,
    centroid,
    area,
    status,
    rng
) where {FT <: AbstractFloat}
    count = 1
    err = FT(1)
    mc_x = zeros(FT, point_generator.npoints)
    mc_y = zeros(FT, point_generator.npoints)
    mc_in = fill(false, point_generator.npoints)
    # Find bounding box
    poly = _translate_poly(FT, poly, -GI.x(centroid), -GI.y(centroid))::Polys{FT}
    (xmin, xmax), (ymin, ymax) = GI.extent(poly)
    Δx, Δy = xmax - xmin, ymax - ymin
    while err > point_generator.err
        if count > 10
            err = 0.0
            status.tag = remove
        else
            mc_x .= xmin .+ Δx * rand(rng, FT, point_generator.npoints)
            mc_y .= ymin .+ Δy * rand(rng, FT, point_generator.npoints)
            mc_in .= [GO.coveredby((mc_x[i], mc_y[i]), poly) for i in eachindex(mc_x)]
            err = abs(sum(mc_in)/point_generator.npoints * (Δx * Δy) - area)/area
            count += 1
        end
    end
    mc_x = mc_x[mc_in]
    mc_y = mc_y[mc_in]
    if isempty(mc_x)
        status.tag = remove
    end

    return mc_x, mc_y, status
end

"""
    generate_subfloe_points(
        point_generator,
        poly,
        centroid,
        area,
        status,
        rng
    )

Generate evenly spaced points within given floe coordinates to be used for
coupling. If only one point falls within the floe, return the floe's centroid.
Inputs:
    point_generator     <SubGridPointsGenerator> sub-grid point generator
    poly                <Polygon> Polygon representing floe shape
    centroid            <Matrix> floe's centroid
    area                <AbstractFloat> floe's area
    status              <Status> floe status (i.e. active, fuse in simulation)
    rng                 <AbstractRNG> random number generator is not used in
                            this generation method
Ouputs:
    x_sub_floe  <Vector{FT}> vector of sub-floe grid points x-coords within floe
    y_sub_floe  <Vector{FT}> vector of sub-floe grid points y-coords within floe
    status      <Status> tag isn't changed with this generation method
"""
function generate_subfloe_points(
    point_generator::SubGridPointsGenerator{FT},
    poly,
    centroid,
    area,
    status,
    rng
) where {FT <: AbstractFloat}
    poly = _translate_poly(FT, poly, -GI.x(centroid), -GI.y(centroid))::Polys{FT}
    (xmin, xmax), (ymin, ymax) = GI.extent(poly)
    xpoints = Vector{FT}()
    ypoints = Vector{FT}()
    local x1, y1
    for (i, p) in enumerate(GI.getpoint(GI.getexterior(poly)))
        # Determine points on edges
        x2, y2 = GI.x(p), GI.y(p)
        if i == 1
            x1, y1 = x2, y2
            continue
        end
        Δx, Δy = x2 - x1, y2 - y1
        l = sqrt((Δx)^2 + (Δy)^2)
        # Add current vertex
        push!(xpoints, x1)
        push!(ypoints, y1)
        # If distance between i and i+1 vertex is less than 2 sub-grid cells
        # but greater than one, add a point inbetween those two vertices
        x2_unshifted, y2_unshifted = x2, y2
        if l <= 2point_generator.Δg
            if l > point_generator.Δg
                push!(xpoints, x1 + Δx/2)
                push!(ypoints, y1 + Δy/2)
            end
        else  # The edge needs more points than the corners and midpoint
            if Δx == 0
                y1 +=  point_generator.Δg/2 * sign(Δy)
                y2 -= point_generator.Δg/2 * sign(Δy)
            elseif Δy == 0
                x1 += point_generator.Δg/2 * sign(Δx)
                x2 -= point_generator.Δg/2 * sign(Δx)
            else  # shift points to still be on the line
                m = Δy / Δx
                x_shift = sqrt(point_generator.Δg^2 / 4(1 + m^2))
                y_shift = m * x_shift
                x1 += x_shift
                x2 -= x_shift
                y1 += y_shift
                y2 -= y_shift
            end
            l = sqrt((x2 - x1)^2 + (y2 - y1)^2)
            n_edge_points = ceil(Int, l / point_generator.Δg) + 1
            append!(xpoints, range(x1, x2, length = n_edge_points))
            append!(ypoints, range(y1, y2, length = n_edge_points))
        end
        x1, y1 = x2_unshifted, y2_unshifted
    end
    # Add points in the interior of the floe
    n_xpoints = ceil(Int, (xmax - xmin) / point_generator.Δg)
    n_ypoints = ceil(Int, (ymax - ymin) / point_generator.Δg)
    x_interior_points = if n_xpoints < 3
        n_xpoints = 1
        FT(0):FT(0)  # polygon is centered at the origin
    else
        range(
            xmin + point_generator.Δg/2,
            xmax - point_generator.Δg/2,
            length = n_xpoints,
        )
    end
    y_interior_points = if n_ypoints < 3
        n_ypoints = 1
        FT(0):FT(0)
    else
        range(
            ymin + point_generator.Δg/2,
            ymax - point_generator.Δg/2,
            length = n_ypoints,
        )
    end
    x_sub_floe = repeat(x_interior_points, n_ypoints)
    y_sub_floe = repeat(y_interior_points, inner = n_xpoints)
    # Check which points are within the polygon and add to list
    in_floe = [GO.coveredby((x_sub_floe[i], y_sub_floe[i]), poly) for i in eachindex(x_sub_floe)]
    append!(xpoints, x_sub_floe[in_floe])
    append!(ypoints, y_sub_floe[in_floe])
    return xpoints, ypoints, status
end

#-------------- Monte Carlo Point Calculations --------------#

"""
    grid_cell_index(p, Δg, g0)

Find the index of given point's cartesian value (in either the x or y direction)
within the simulation grid. 
Inputs:
    p       <Real> point's cartesian value in either x or y direction
    Δg      <Real> simulation grid's cell width or height
    g0      <Real> simulation grid's first grid line value in either x or y
                        direction
Output:
    Point's grid cell index within the simulation grid, as specified by the grid
    cell dimension and grid line starting value, in either the x or y direction
"""
grid_cell_index(p, Δg, g0) = floor(Int, (p - g0)/Δg) + 1

"""
    grid_line_index(p, Δg, g0)

Find the index of given point's cartesian value (in either the x or y direction)
within grid with cells centered on simulation grid's grid lines. Thus these
cells are shifted from simulation's grid cells by half of a grid cell to the
left.  
Inputs:
    p       <Real> point's cartesian value in either x or y direction
    Δg      <Real> grid's cell width or height
    g0      <Real> grid's first grid line value in either x or y direction
Output:
    Point's grid cell index within the shifted simulation grid, as specified by
    the grid cell dimension and grid line starting value, in either the x or y
    direction
"""
grid_line_index(p, Δg, g0) = floor(Int, (p - g0)/Δg + 0.5) + 1

"""
    grid_xg_index(xp, yp, grid::RegRectilinearGrid)

Find indices of given cartesian point within simulation's xg-grid.
Inputs:
    xp      <Real> point's x-cartesian value
    yp      <Real> point's y-cartesian value
    grid    <RegRectilinearGrid> simulation's grid
Outputs:
    x index and y indices within xg-grid
"""
grid_xg_index(xp, yp, grid::RegRectilinearGrid) = 
    grid_line_index(xp, grid.Δx, grid.x0),
    grid_cell_index(yp, grid.Δy, grid.y0)

"""
    grid_yg_index(xp, yp, grid::RegRectilinearGrid)

Find indices of given cartesian point within simulation's yg-grid.
Inputs:
    xp      <Real> point's x-cartesian value
    yp      <Real> point's y-cartesian value
    grid    <RegRectilinearGrid> simulation's grid
Outputs:
    x index and y indices within yg-grid
"""
grid_yg_index(xp, yp, grid::RegRectilinearGrid) = 
    grid_cell_index(xp, grid.Δx, grid.x0),
    grid_line_index(yp, grid.Δy, grid.y0)

"""
    grid_xc_index(xp, yp, grid::RegRectilinearGrid)

Find indices of given cartesian point within simulation's xc-grid.
Inputs:
    xp      <Real> point's x-cartesian value
    yp      <Real> point's y-cartesian value
    grid    <RegRectilinearGrid> simulation's grid
Outputs:
    x index and y indices within xc-grid
Note: 
    This is equivalent fo the yc-grid
"""   
grid_xc_index(xp, yp, grid::RegRectilinearGrid) =
    grid_cell_index(xp, grid.Δx, grid.x0),
    grid_cell_index(yp, grid.Δy, grid.y0)

"""
    grid_yc_index(xp, yp, grid::RegRectilinearGrid)

Find indices of given cartesian point within simulation's yc-grid.
Inputs:
    xp      <Real> point's x-cartesian value
    yp      <Real> point's y-cartesian value
    grid    <RegRectilinearGrid> simulation's grid
Outputs:
    x index and y indices within yc-grid
Note: 
    This is equivalent fo the xc-grid
"""   
grid_yc_index(xp, yp, grid::RegRectilinearGrid) =
    grid_xc_index(xp, yp, grid)

"""
    find_grid_cell_index(xp, yp, grid::RegRectilinearGrid)
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
"""
function find_grid_cell_index(xp, yp, grid::RegRectilinearGrid)
    xidx = floor(Int, (xp - grid.x0) / grid.Δx) + 1
    yidx = floor(Int, (yp - grid.y0) / grid.Δy) + 1
    return xidx, yidx
end


"""
    find_center_cell_index(xp, yp, grid::RegRectilinearGrid)
Find index of the cell centered on grid lines of the given RegRectilinearGrid
that the given x-coordinate and y-coordinate falls within.
This cell is centered around the grid lines, so it is a shifted grid cell
by half a cell. Method depends on grid being a regular rectilinear grid.
Inputs:
    xp      <AbstractFloat> x-coordinates of point
    yp      <AbstractFloat> y-coordinate of point
    grid    <RegRectilinearGrid> simulation grid
Outputs:
    xidx    <AbstractFloat> x-index of grid cell (cented on grid lines) x-point
                is within - this is the column
    yidx    <AbstractFloat> y-index of grid cell (cented on grid lines) y-point
                is within - this is the row
Note:
    Points can be outside of the grid, so index can be less than 1 or greater
    than the number of grid lines in a given direction.
"""
function find_center_cell_index(xp, yp, grid::RegRectilinearGrid)
    xidx = floor(Int, (xp - grid.x0)/(grid.Δx) + 0.5) + 1
    yidx = floor(Int, (yp - grid.y0)/(grid.Δy) + 0.5) + 1
    return xidx, yidx
end
"""
    in_bounds(
        xr,
        yr,
        grid,
        ::NonPeriodicBoundary,
        ::NonPeriodicBoundary,
    )

With all non-periodic boundaries, points outside of the grid in both the x and y
are defined to be out of bounds since these points can't be interpolated as we
don't have any information on the ocean outside of the grid.
Inputs:
    xr  <AbstractFloat> point x-coordinate
    yr  <AbstractFloat> point y-coordinate
        <::NonPeriodicBoundary> type of either north or south boundary -
            checking if periodic pair
        <::NonPeriodicBoundary> type of either east or west boundary -
            checking if periodic pair
Output:
    Boolean that is true if both xr and yr are within domain boundaries, and
    false otherwise.
"""
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

"""
    function in_bounds(
        xr,
        yr,
        grid,
        ::NonPeriodicBoundary,
        ::PeriodicBoundary,
    )

With the north/south non-periodic boundaries, points outside of the grid in the
y-direction are defined to be out of bounds since these points can't be
interpolated as we don't have any information on the ocean outside of the grid.
Inputs:
    xr  <AbstractFloat> point x-coordinate
    yr  <AbstractFloat> point y-coordinate
        <::NonPeriodicBoundary> type of either north or south boundary -
            checking if periodic pair
        <::PeriodicBoundary> type of either east or west boundary -
            checking if periodic pair
Output:
    Boolean that is true if yr is within domain boundaries, and false otherwise.
"""
function in_bounds(
    xr,
    yr,
    grid,
    ::NonPeriodicBoundary,
    ::PeriodicBoundary,
)
    return grid.y0 <= yr <= grid.yf
end

"""
    function in_bounds(
        xr,
        yr,
        grid,
        ::PeriodicBoundary,
        ::NonPeriodicBoundary,
    )

With the east/west non-periodic boundaries, points outside of the grid in the
x-direction are defined to be out of bounds since these points can't be
interpolated as we don't have any information on the ocean outside of the grid.
Inputs:
    xr  <AbstractFloat> point x-coordinate
    yr  <AbstractFloat> point y-coordinate
        <::PeriodicBoundary> type of either north or south boundary -
            checking if periodic pair
        <::NonPeriodicBoundary> type of either east or west boundary -
            checking if periodic pair
Output:
    Boolean that is true if xr is within domain boundaries, and false otherwise.
"""
function in_bounds(
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::NonPeriodicBoundary,
)
    return grid.x0 <= xr <= grid.xf
end

"""
    function in_bounds(
        xr,
        yr,
        grid,
        ::PeriodicBoundary,
        ::PeriodicBoundary,
    )

With all periodic boundaries, all points are considered to be in-bounds.
Inputs:
    xr  <AbstractFloat> point x-coordinate
    yr  <AbstractFloat> point y-coordinate
        <::PeriodicBoundary> type of either north or south boundary -
            checking if periodic pair
        <::PeriodicBoundary> type of either east or west boundary -
            checking if periodic pair
Output:
    Boolean that is true regardless of point values.
"""
function in_bounds(
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::PeriodicBoundary,
)
    return true
end

"""
    calc_subfloe_values!(
        floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
        grid,
        domain,
        mc_cart,
        mc_grid_idx,
    )

Calculates subfloe point's cartesian coordiantes, polar coordiantes,
velocity and index within the grid. 
Inputs:
    floe        <Union{Floe{AbstractFloat}, LazyRow{Floe{AbstractFloat}}}> floe
    grid        <AbstractRectilinearGrid> model's grid
    domain      <Domain> model's domain
    mc_cart     <Matrix{AbstractFloat}> pre-allocated nx2 matrix for floe's
                    monte carlo point's cartesian coordinates where the first
                    column is x and second is y
    mc_grid_idx <Matrix{AbstractFloat}> pre-allocated nx2 matrix for floe's
                    monte carlo point's grid indices where the first column is
                    the column and the second is the row that the point is in on
                    the grid split into cells centered on grid lines.
Outputs:
    j   <Int> last element in mc_cart and mc_grid_idx that holds monte carlo
            point information for given floe.
    mc_cart and mc_grid_idx filled with data for given floe's monte carlo points
    up to row j.
"""
function calc_subfloe_values!(
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
        if in_bounds(x, y, grid, domain.north, domain.east)
            j += 1  # if added to outputs, move to next index in output array
            cart_vals[j, 1] = x
            cart_vals[j, 2] = y
            grid_idx[j, 1], grid_idx[j, 2] = find_center_cell_index(
                cart_vals[j, 1],
                cart_vals[j, 2],
                grid,
            )
        end
    end
    return j  # last spot in output arrays that has values for this floe
end

#-------------- Interpolation of Ocean and Atmosphere --------------#

"""
    find_interp_knots(
        point_idx,
        ncells,
        L,
        Δd::Int,
        ::PeriodicBoundary,
    )

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
"""
function find_interp_knots(
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

"""
    find_interp_knots(
        point_idx,
        ncells,
        L,
        Δd::Int,
        ::NonPeriodicBoundary,
    )

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
"""
function find_interp_knots(
    point_idx,
    ncells,
    glines,
    L,
    Δd::Int,
    ::NonPeriodicBoundary,
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

"""
    mc_interpolation(
        mc_cart,
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
    mc_cart             <Matrix{AbstractFloat}> cartesian coordinates for model
                            coordinates - nx2 matrix of monte carlo coordinates
                            where first column is the x-coords and the second
                            column is the y-coords
    mc_grid_idx         <Matrix{Int}> index of monte carlo points within the
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
"""
function mc_interpolation(
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
    xknots, xknot_idx = find_interp_knots(
        xidx,
        grid.Nx,
        grid.x0:grid.Δx:grid.xf,
        grid.xf - grid.x0,
        coupling_settings.Δd,
        domain.east,
    )
    yknots, yknot_idx = find_interp_knots(
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

#-------------- Effects of Ice and Atmosphere on Ocean --------------#
"""
    check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ::PeriodicBoundary,
        ::PeriodicBoundary,
    )

Return cell bounding values as is given the domain is doubley periodic and thus
the cell can extend beyond the grid as it will simply wrap back around into grid
through opposite periodic boundary.
Inputs:
    xmin    <Float> center cell minimum x value
    xmax    <Float> center cell maxumum x value
    ymin    <Float> center cell minimum y value
    ymax    <Float> center cell maximum y value
    grid    <AbstractRectilinearGrid> model's grid
            <PeriodicBoundary> type of north or south boundary - periodic pair
            <PeriodicBoundary> type of east or west boundary - periodic pair
Output:
    x and y minimums and maximums as given since they can extend past the grid
    due to periodic boundaries.
"""
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

"""
    check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ::NonPeriodicBoundary,
        ::PeriodicBoundary,
    )

Trim cell bound in the north-south direction if it exends past grid due to
non-periodic boundary pair.
Inputs:
    xmin    <Float> center cell minimum x value
    xmax    <Float> center cell maxumum x value
    ymin    <Float> center cell minimum y value
    ymax    <Float> center cell maximum y value
    grid    <AbstractRectilinearGrid>
            <NonPeriodicBoundary> type of either north or south boundary - not a
                periodic pair
            <PeriodicBoundary> type of either east or west boundary - periodic
                pair
Output:
        Potentially trimmed y min and y max if these values extend beyond grid
        values. Else returned unchanged. 
"""
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

"""
    check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ::PeriodicBoundary,
        ::NonPeriodicBoundary,
    )

Trim cell bound in the east-west direction if it exends past grid due to
non-periodic boundary pair.
Inputs:
    xmin    <Float> center cell minimum x value
    xmax    <Float> center cell maxumum x value
    ymin    <Float> center cell minimum y value
    ymax    <Float> center cell maximum y value
    grid    <AbstractRectilinearGrid>
            <PeriodicBoundary> type of either north or south boundary - periodic
                pair
            <NonPeriodicBoundary> type of either east or west boundary - not a
                periodic pair
Output:
    Potentially trimmed x min and max if these values extend beyond grid values.
    Else returned unchanged. 
"""
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

"""
    check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ::NonPeriodicBoundary,
        ::NonPeriodicBoundary,
    )

Trim cell bounds in the east-west and north-south direction if they exend past
grid due to non-periodic boundary pairs.
Inputs:
    xmin    <Float> center cell minimum x value
    xmax    <Float> center cell maxumum x value
    ymin    <Float> center cell minimum y value
    ymax    <Float> center cell maximum y value
    grid    <AbstractRectilinearGrid>
            <NonPeriodicBoundary> type of either north or south boundary - not a
                periodic pair
            <NonPeriodicBoundary> type of either east or west boundary - not a
                periodic pair
Output:
    Potentially trimmed x and y minimums and maximums if these values extend
    beyond grid values. Else returned unchanged. 
"""
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

"""
    center_cell_coords(
        xidx::Int,
        yidx::Int,
        grid::RegRectilinearGrid,
        ns_bound,
        ew_bound,
    )

Find the coordinates of a given grid cell, centered on a grid line with row yidx
and column xidx. This is offset from the cells within the regular rectilinear
grid by half of a grid cell. 
Inputs:
    xidx        <Int> x index of grid line within list of gridlines (cell column)
    yidx        <Int> y index of grid line within list of gridlines (cell row)
    grid        <RegRectilinearGrid> model's grid 
    ns_bound    <AbstractBoundary> type of either north or south boundary - for
                    checking if periodic
    ew_bound    <AbstractBoundary> type of either east or west boundary - for
                    checking if perioidic
Output:
    <PolyVec> coordinates for cell centered on grid line with given indices.
    Note that cell bounds will be adjusted depending on if the bounds are
    periodic. Cells cannot extend outside of non-periodic boundaries and thus
    will be trimmed at boundaries. Therefore, if indices place cell completely
    outside of grid, could return a line at the edge of the boundary. 
"""
function center_cell_coords(
    ::Type{FT},
    xidx::Int,
    yidx::Int,
    grid::RegRectilinearGrid,
    ns_bound,
    ew_bound,
) where FT
    xmin = (xidx - 1.5) * grid.Δx + grid.x0
    xmax = xmin + grid.Δx
    ymin = (yidx - 1.5) * grid.Δy + grid.y0
    ymax = ymin + grid.Δy
    #= Check if cell extends beyond boundaries and if non-periodic, trim cell to
    fit within grid. =#
    xmin, xmax, ymin, ymax = check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ns_bound,
        ew_bound,
    )
    return _make_bounding_box_polygon(FT, xmin, xmax, ymin, ymax)
end

"""
    shift_cell_idx(idx, nlines, ::NonPeriodicBoundary)

Return index as is given non-periodic boundary pair in either x or y direction.
Inputs:
    idx     <Int> grid line index in either x or y
    nlines  <Int> number of grid lines in model grid in either x or y direction
            <NonPeriodicBoundary> boundary pair is non-periodic
Ouput:
    idx <Int> as given. Can include the index nlines, unlike with the periodic
        case, which will use the first index instead. 
"""
function shift_cell_idx(idx, nlines, ::NonPeriodicBoundary)
    return idx
end

"""
    shift_cell_idx(idx, nlines, ::PeriodicBoundary)

If index is greater than or equal to the grid lines, shift index to equivalent
grid line on opposite side of grid due to periodic boundary. Similarly if given
index is less than 1, shift index to equivalent grid line on opposite side of
grid due to periodic boundary.
Inputs:
    idx     <Int> grid line index in either x or y
    nlines  <Int> number of grid lines in model grid in either x or y direction
            <PeriodicBoundary> boundary pair is periodic
Output:
    <Int> if given index is greater than or equal to number of grid lines, shift
    index. If given index is less than 1, shift grid index. For example, the
    last grid index, nlines, is equivalent to the 1st grid line. The nlines+1
    grid line is equivalent to the 2nd grid line.
"""
function shift_cell_idx(idx, nlines, ::PeriodicBoundary)
    ncells = nlines - 1
    return idx < 1 ? (idx + ncells) : ncells < idx ? (idx - ncells) : idx
end

#-------------- Ocean and Atmosphere on Ice --------------#
"""
    calc_atmosphere_forcing(
        mc_xr, 
        mc_yr,
        upoint,
        vpoint,
        uatm_interp,
        vatm_interp,
        c,
    )

Calculates the stresses on a floe from the atmosphere above at given monte
carlo point.
Inputs:
    mc_xr       <AbstractFloat> monte carlo point x-coordinate
    mc_yr       <AbstractFloat> monte carlo point y-coordinate
    upoint      <AbstractFloat> u velocity of floe at monte carlo point
    vpoint      <AbstractFloat> v velocity of floe at monte carlo point
    uatm_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the atompshere u velocity onto point
    vatm_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the atompshere v velocity onto point
    c           <Constants> simulation's constants
Outputs:
    τx_atm  <AbstractFloat> stress from atmosphere on floe in
                x-direction at given monte carlo point
    τy_atm  <AbstractFloat> stress from atmosphere on floe in
                y-direction at given monte carlo point
"""
function calc_atmosphere_forcing(
    xr, 
    yr,
    upoint,
    vpoint,
    uatm_interp,
    vatm_interp,
    c,  # constants
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

"""
    calc_ocean_forcing!(
        mc_xr,
        mc_yr,
        upoint,
        vpoint,
        uocn_interp,
        vocn_interp,
        hflx_interp,
        ma_ratio,
        c,
    )

Calculates the stresses on a floe from the ocean above at given monte carlo
point.
Inputs:
    mc_xr       <AbstractFloat> monte carlo point x-coordinate
    mc_yr       <AbstractFloat> monte carlo point y-coordinate
    upoint      <AbstractFloat> u velocity of floe at monte carlo point
    vpoint      <AbstractFloat> v velocity of floe at monte carlo point
    uocn_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean u velocity onto point
    vocn_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean v velocity onto point
    hflx_interp <Interplations object> linear interpolation function from
                    Interpolations.jl that takes in two arguments (x, y) and
                    interpolates the ocean heatflux factor velocity onto point
    ma_ratio    <AbstractFloat> floe's mass to area ratio
    c           <Constants> simulation's constants
Outputs:
    τx_ocn          <AbstractFloat> stress from ocean velocity on floe in
                        x-direction at given monte carlo point
    τy_ocn          <AbstractFloat> stress from ocean velocity on floe in
                        y-direction at given monte carlo point
    τx_pressure∇    <AbstractFloat> stress from ocean pressure gradient on floe
                        in x-direction at given monte carlo point
    τy_pressure∇    <AbstractFloat> stress from ocean pressure gradient on floe
                        in y-direction at given monte carlo point
    hflx_factor     <AbstractFloat> heatflux factor at given monte carlo point
                        from the heatflux factors of ocean below floe
"""
function calc_ocean_forcing!(
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

"""
    add_point!(
        cfloes::CellFloes,
        scell::CellStresses,
        floeidx,
        τx,
        τy,
        Δx,
        Δy,
    )
Add floe to CellFloes list of floes within that grid cell and aggragate the
stress caused by monte carlo point in floe into CellStresses object.
Inputs:
    cfloes  <CellFloes> CellFloes object representing one grid cell (centered on
                model's grid lines)
    scell   <CellStresses> CellStresses aggragating stresses from floes within
                grid cell from each floes' monte carlo points
    floeidx <Int> floe index within model's list of floes
    τx      <AbstractFloat> x-directional stress from monte carlo point on ocean
    τy      <AbstractFloat> y-directional stress from monte carlo point on ocean
    Δx      <AbstractFloat> x-translation to move floe from current position into
                given grid cell if shifted due to periodic boundaries
    Δy      <AbstractFloat> y-translation to move floe from current position into
                given grid cell if shifted due to periodic boundaries
Outputs:
    None. Add information to both cfloes and scell to aggregate stress on ocean
    grid cell and record where floe is on model grid. 
"""
function add_point!(
    cfloes::CellFloes,
    scell::CellStresses,
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
"""
    add_point!(
        cfloes::CellFloes,
        floeidx,
        Δx,
        Δy,
    )
Add floe to CellFloes list of floes within that grid cell and aggragate the
stress caused by monte carlo point in floe into CellStresses object.
Inputs:
    cfloes  <CellFloes> CellFloes object representing one grid cell (centered on
                model's grid lines)
    floeidx <Int> floe index within model's list of floes
    Δx      <AbstractFloat> x-translation to move floe from current position into
                given grid cell if shifted due to periodic boundaries
    Δy      <AbstractFloat> y-translation to move floe from current position into
                given grid cell if shifted due to periodic boundaries
Outputs:
    None. Add information to both cfloes to record where floe is on model grid. 
"""
function add_point!(
    cfloes::CellFloes,
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

"""
    floe_to_grid_info!(
        floeidx,
        xidx,
        yidx,
        τx_ocn::FT,
        τy_ocn::FT,
        grid,
        domain,
        scells,
    )

Add force from the ice on ocean to ocean force fields (fx & fy) for each grid
cell and update ocean sea ice area fraction (si_area), representing total area
of sea ice in a given cell. Function is called for each monte carlo point.
Inputs:
    floeidx             <Int> index of floe within model's floe array
    xidx                <Int> grid x index that floe's point is within for grid
                            centered on grid lines
    yidx                <Int> grid column that floe's point is within for grid
                            centered on grid lines
    τx_ocn              <AbstractFloat> x-stress caused by ocean on point
    τy_ocn              <AbstractFloat> y-stress caused by ocean on point
    grid                <AbstractRectilinearGrid> model's grid
    domain              <Domain> model's domain
    cell_floes          <Matrix{CellFloes}> matrix of CellFloes, one for each
                            grid cell
    scells              <Matrix{CellStresses}> matrix of CellStressess, one
                            for each grid cell
    coupling_settings   <CouplingSettings> simulation's coupling settings
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
    scells,
    coupling_settings,
) where {FT}
    # Determine grid cell point is in and if floe is shifted by periodic bounds
    shifted_xidx = shift_cell_idx(xidx, grid.Nx + 1, ew_bound)
    shifted_yidx = shift_cell_idx(yidx, grid.Ny + 1, ns_bound)
    Δx = (shifted_xidx - xidx) * (grid.Δx)
    Δy = (shifted_yidx - yidx) * (grid.Δy)
    if coupling_settings.two_way_coupling_on 
        # If two-way coupling, save stress on ocean per cell
        add_point!(
            grid.floe_locations[shifted_xidx, shifted_yidx],
            scells[shifted_xidx, shifted_yidx],
            floeidx,
            -τx_ocn,
            -τy_ocn,
            Δx,
            Δy,
        )
    else
        add_point!(
            grid.floe_locations[shifted_xidx, shifted_yidx],
            floeidx,
            Δx,
            Δy,
        )
    end
    return
end

"""
    calc_one_way_coupling!(
        floes::StructArray{Floe{FT}},
        grid,
        atmos,
        ocean,
        domain,
        coupling_settings,
        consts,
    )

Preforms calculations needed for one way coupling by calculating floe's forcings
from ocean and atmosphere as well as the heatflux below a given floe.

Floe location on grid is also recorded. If two-way coupling is on, total
stress on each grid cell per-floe in grid cell is also recorded for use in
calc_two_way_coupling!
Inputs:
    floes               <StructArray{Floe{FT}}> model's floe list
    grid                <AbstractRectilinearGrid> model's grid
    atmos               <Ocean> model's atmosphere
    ocean               <Ocean> model's ocean
    domain              <Domain> model's domain
    coupling_settings   <CouplingSettings> simulation coupling settings
    consts              <Constants> simulation's constants
Ouputs:
    None. Update each floe's forces, torque, and heatflux factor from
    ocean/atmosphere. Determine location of floe within grid and if two-way
    coupling in enabled, save floe stress on grid. 
"""
function calc_one_way_coupling!(
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
        npoints = calc_subfloe_values!(
            get_floe(floes, i),
            grid,
            domain,
            cart_vals,
            grid_idx,
        )
        if npoints == 0
            floes.status[i].tag = remove
        else
            # Interpolaters for ocean and atmosphere
            uatm_int, vatm_int, uocn_int, vocn_int, hflx_int = mc_interpolation(
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
                τx_atm, τy_atm = calc_atmosphere_forcing(
                    cart_vals[j, 1], 
                    cart_vals[j, 2],
                    upoint,
                    vpoint,
                    uatm_int,
                    vatm_int,
                    consts,
                )
                τx_ocn, τy_ocn, τx_p∇, τy_p∇, hflx_factor = calc_ocean_forcing!(
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
                floe_to_grid_info!(
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

"""
    calc_two_way_coupling!(
        floes::StructArray{Floe{FT}},
        grid::RegRectilinearGrid,
        atmos,
        ocean,
        domain,
        floe_settings,
        consts,
        Δt,
    )

Calculate effects of ice and atmosphere on the ocean and update ocean stress
fields and sea ice fraction.
Inputs:
    floes           <StructArray{Floe}> model's floes
    grid            <AbstractRectilinearGrid> model's grid
    atmos           <Atmos> model's atmosphere
    ocean           <Ocean> model's ocean
    domain          <Domain> model's domain
    floe_settings   <FloeSettings> simulation's floe settings
    consts          <Constants> model's constants
    Δt              <Int> simulation's timestep in seconds
Output:
    None. Update's ocean's stress fields and heatflux factor field. 
"""
function calc_two_way_coupling!(
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
    Threads.@threads  for cartidx in CartesianIndices(ocean.scells)
        ocean.τx[cartidx] = FT(0)
        ocean.τy[cartidx] = FT(0)
        ocean.si_frac[cartidx] = FT(0)
        τocn = ocean.scells[cartidx]
        floe_locations = grid.floe_locations[cartidx]
        if !isempty(floe_locations.floeidx)
            # Coordinates of grid cell
            cell_poly = center_cell_coords(
                FT,
                cartidx[1],
                cartidx[2],
                grid,
                domain.north,
                domain.east
            )
            for i in eachindex(floe_locations.floeidx)
                floe_poly = _translate_poly(FT,
                    floes.poly[floe_locations.floeidx[i]],
                    floe_locations.Δx[i],
                    floe_locations.Δy[i],
                )::Polys{FT}
                floe_area_in_cell = sum(
                    GO.area.(intersect_polys(cell_poly, floe_poly), FT)
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


"""
    timestep_coupling!(
        model,
        Δt,
        consts,
        coupling_settings,
        floe_settings,
    )

Calculates the effects of the ocean and atmosphere on the ice and the effects of
the ice and atmosphere on the ocean if the coupling is two-way.
Inputs:
    model               <Model> model
    Δt                  <Int> length of timestep in seconds
    consts              <Constants> constants used in simulation
    coupling_settings   <CouplingSettings> settings for coupling
    floe_settings       <FloeSettings> settings for basic floe properties
Outputs:
    None. Updates each floe's ocean/atmosphere forcings (fxOA, fyOA, torqueOA)
    and calculates stresses on each ocean grid cell from ice and atmosphere if
    two-way coupling is enabled in coupling_settings       
"""
function timestep_coupling!(
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
    calc_one_way_coupling!(
        model.floes,
        model.grid,
        model.atmos,
        model.ocean,
        model.domain,
        coupling_settings,
        consts,
    )
    if coupling_settings.two_way_coupling_on
        calc_two_way_coupling!(
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
